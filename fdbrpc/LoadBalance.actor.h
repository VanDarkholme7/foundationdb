/*
 * LoadBalance.actor.h
 *
 * This source file is part of the FoundationDB open source project
 *
 * Copyright 2013-2018 Apple Inc. and the FoundationDB project authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

// When actually compiled (NO_INTELLISENSE), include the generated version of this file.  In intellisense use the source
// version.
#if defined(NO_INTELLISENSE) && !defined(FLOW_LOADBALANCE_ACTOR_G_H)
#define FLOW_LOADBALANCE_ACTOR_G_H
#include "fdbrpc/LoadBalance.actor.g.h"
#elif !defined(FLOW_LOADBALANCE_ACTOR_H)
#define FLOW_LOADBALANCE_ACTOR_H

#include "flow/flow.h"
#include "flow/Knobs.h"

#include "fdbrpc/FailureMonitor.h"
#include "fdbrpc/fdbrpc.h"
#include "fdbrpc/Locality.h"
#include "fdbrpc/QueueModel.h"
#include "fdbrpc/MultiInterface.h"
#include "flow/actorcompiler.h" // This must be the last #include.

using std::vector;

struct ModelHolder : NonCopyable, public ReferenceCounted<ModelHolder> {
	QueueModel* model;
	bool released;
	double startTime;
	double delta;
	uint64_t token;

	ModelHolder(QueueModel* model, uint64_t token) : model(model), token(token), released(false), startTime(now()) {
		if (model) {
			delta = model->addRequest(token);
		}
	}

	void release(bool clean, bool futureVersion, double penalty, bool measureLatency = true) {
		if (model && !released) {
			released = true;
			double latency = (clean || measureLatency) ? now() - startTime : 0.0;
			model->endRequest(token, latency, penalty, delta, clean, futureVersion);
		}
	}

	~ModelHolder() { release(false, false, -1.0, false); }
};

// Subclasses must initialize all members in their default constructors
// Subclasses must serialize all members
struct LoadBalancedReply {
	double penalty;
	Optional<Error> error;
	LoadBalancedReply() : penalty(1.0) {}
};

Optional<LoadBalancedReply> getLoadBalancedReply(const LoadBalancedReply* reply);
Optional<LoadBalancedReply> getLoadBalancedReply(const void*);

// Stores state for a request made by the load balancer
template <class Request>
struct RequestData : NonCopyable {
	typedef ErrorOr<REPLY_TYPE(Request)> Reply;

	Future<Reply> response;
	Reference<ModelHolder> modelHolder;
	bool triedAllOptions = false;

	bool requestStarted = false; // true once the request has been sent to an alternative
	bool requestProcessed = false; // true once a response has been received and handled by checkAndProcessResult

	// Whether or not the response future is valid
	// This is true once setupRequest is called, even though at that point the response is Never().
	bool isValid() { return response.isValid(); }

	// Initializes the request state and starts it, possibly after a backoff delay
	void startRequest(double backoff,//延迟时间
	                  bool triedAllOptions,//是否尝试了所有选项
	                  RequestStream<Request> const* stream,//请求流
	                  Request const& request,//具体的请求
	                  QueueModel* model) {
		modelHolder = Reference<ModelHolder>();
		requestStarted = false;

		if (backoff > 0) {
			response = mapAsync<Void, std::function<Future<Reply>(Void)>, Reply>(
			    delay(backoff), [this, stream, &request, model](Void _) {
				    requestStarted = true;
				    modelHolder = Reference<ModelHolder>(new ModelHolder(model, stream->getEndpoint().token.first()));
				    return stream->tryGetReply(request);
			    });
		} else {
			requestStarted = true;
			modelHolder = Reference<ModelHolder>(new ModelHolder(model, stream->getEndpoint().token.first()));
			//发送请求，等待响应
			response = stream->tryGetReply(request);
		}
		//更新请求状态
		requestProcessed = false;
		this->triedAllOptions = triedAllOptions;
	}

	// Implementation of the logic to handle a response.
	// Checks the state of the response, updates the queue model, and returns one of the following outcomes:
	// A return value of true means that the request completed successfully
	// A return value of false means that the request failed but should be retried
	// A return value with an error means that the error should be thrown back to original caller
	static ErrorOr<bool> checkAndProcessResultImpl(Reply const& result,
	                                               Reference<ModelHolder> modelHolder,
	                                               bool atMostOnce,
	                                               bool triedAllOptions) {
		ASSERT(modelHolder);

		Optional<LoadBalancedReply> loadBalancedReply;
		if (!result.isError()) {
			loadBalancedReply = getLoadBalancedReply(&result.get());
		}

		int errCode;
		if (loadBalancedReply.present()) {
			errCode = loadBalancedReply.get().error.present() ? loadBalancedReply.get().error.get().code()
			                                                  : error_code_success;
		} else {
			errCode = result.isError() ? result.getError().code() : error_code_success;
		}

		bool maybeDelivered = errCode == error_code_broken_promise || errCode == error_code_request_maybe_delivered;
		bool receivedResponse =
		    loadBalancedReply.present() ? !loadBalancedReply.get().error.present() : result.present();
		receivedResponse = receivedResponse || (!maybeDelivered && errCode != error_code_process_behind);
		bool futureVersion = errCode == error_code_future_version || errCode == error_code_process_behind;

		modelHolder->release(
		    receivedResponse, futureVersion, loadBalancedReply.present() ? loadBalancedReply.get().penalty : -1.0);

		if (errCode == error_code_server_overloaded) {
			return false;
		}

		if (loadBalancedReply.present() && !loadBalancedReply.get().error.present()) {
			return true;
		}

		if (!loadBalancedReply.present() && result.present()) {
			return true;
		}

		if (receivedResponse) {
			return loadBalancedReply.present() ? loadBalancedReply.get().error.get() : result.getError();
		}

		if (atMostOnce && maybeDelivered) {
			return request_maybe_delivered();
		}

		if (triedAllOptions && errCode == error_code_process_behind) {
			return process_behind();
		}

		return false;
	}

	// Checks the state of the response, updates the queue model, and returns one of the following outcomes:
	// A return value of true means that the request completed successfully
	// A return value of false means that the request failed but should be retried
	// In the event of a non-retryable failure, an error is thrown indicating the failure
	bool checkAndProcessResult(bool atMostOnce) {
		ASSERT(response.isReady());
		requestProcessed = true;

		ErrorOr<bool> outcome =
		    checkAndProcessResultImpl(response.get(), std::move(modelHolder), atMostOnce, triedAllOptions);

		if (outcome.isError()) {
			throw outcome.getError();
		} else if (!outcome.get()) {
			response = Future<Reply>();
		}

		return outcome.get();
	}

	// Convert this request to a lagging request. Such a request is no longer being waited on, but it still needs to be
	// processed so we can update the queue model.
	void makeLaggingRequest() {
		ASSERT(response.isValid());
		ASSERT(!response.isReady());
		ASSERT(modelHolder);
		ASSERT(modelHolder->model);

		QueueModel* model = modelHolder->model;
		if (model->laggingRequestCount > FLOW_KNOBS->MAX_LAGGING_REQUESTS_OUTSTANDING ||
		    model->laggingRequests.isReady()) {
			model->laggingRequests.cancel();
			model->laggingRequestCount = 0;
			model->addActor = PromiseStream<Future<Void>>();
			model->laggingRequests = actorCollection(model->addActor.getFuture(), &model->laggingRequestCount);
		}

		// We need to process the lagging request in order to update the queue model
		Reference<ModelHolder> holderCapture = std::move(modelHolder);
		bool triedAllOptionsCapture = triedAllOptions;
		Future<Void> updateModel = map(response, [holderCapture, triedAllOptionsCapture](Reply result) {
			checkAndProcessResultImpl(result, holderCapture, false, triedAllOptionsCapture);
			return Void();
		});
		model->addActor.send(updateModel);
	}

	~RequestData() {
		// If the request has been started but hasn't completed, mark it as a lagging request
		if (requestStarted && !requestProcessed && modelHolder && modelHolder->model) {
			makeLaggingRequest();
		}
	}
};

// Try to get a reply from one of the alternatives until success, cancellation, or certain errors.
// Load balancing has a budget to race requests to a second alternative if the first request is slow.
// Tries to take into account failMon's information for load balancing and avoiding failed servers.
// If ALL the servers are failed and the list of servers is not fresh, throws an exception to let the caller refresh the
// list of servers.
// When model is set, load balance among alternatives in the same DC aims to balance request queue length on these
// interfaces. If too many interfaces in the same DC are bad, try remote interfaces.
ACTOR template <class Interface, class Request, class Multi>
Future<REPLY_TYPE(Request)> loadBalance(
    Reference<MultiInterface<Multi>> alternatives,//备选服务器集合
    RequestStream<Request> Interface::*channel, //请求的通信通道
    Request request = Request(), //要发送的实际请求
    TaskPriority taskID = TaskPriority::DefaultPromiseEndpoint,//任务优先级
    bool atMostOnce = false, // if true, throws request_maybe_delivered() instead of retrying automatically 请求只发送一次
    QueueModel* model = nullptr) { //队列模型，默认为空，用于负载均衡

	//状态变量，用于管理第一个和第二个请求的数据结构
	state RequestData<Request> firstRequestData;
	state RequestData<Request> secondRequestData;

	//第一个请求的终端点标识符
	state Optional<uint64_t> firstRequestEndpoint;
	state Future<Void> secondDelay = Never();

	state Promise<Void> requestFinished;//标记请求是否完成
	state double startTime = now();
	//设置请求优先级
	setReplyPriority(request, taskID);

	//确保至少有1个可选服务器
	if (!alternatives)
		return Never();
	ASSERT(alternatives->size());

	//随机选择初始的最佳服务器和next备选服务器
	state int bestAlt = deterministicRandom()->randomInt(0, alternatives->countBest());
	state int nextAlt = deterministicRandom()->randomInt(0, std::max(alternatives->size() - 1, 1));
	//确保这两个服务器不重合
	if (nextAlt >= bestAlt)
		nextAlt++;
	//如果提供了队列模型，计算每个备选服务器的负载均衡指标
	//选择最佳和次佳的备选服务器，并确定是否需要延迟向次佳服务器发送请求
	if (model) {
		double bestMetric = 1e9; // Storage server with the least outstanding requests.
		double nextMetric = 1e9;
		double bestTime = 1e9; // The latency to the server with the least outstanding requests.
		double nextTime = 1e9;
		int badServers = 0;

		//遍历所有服务器，根据请求队列长度、延迟等指标选择最佳和次佳服务器
		for (int i = 0; i < alternatives->size(); i++) {
			// countBest(): the number of alternatives in the same locality (i.e., DC by default) as alternatives[0].
			// if the if-statement is correct, it won't try to send requests to the remote ones.
			if (badServers < std::min(i, FLOW_KNOBS->LOAD_BALANCE_MAX_BAD_OPTIONS + 1) &&
			    i == alternatives->countBest()) {
				// When we have at least one healthy local server, and the bad
				// server count is within "LOAD_BALANCE_MAX_BAD_OPTIONS". We
				// do not need to consider any remote servers.
				break;
			}

			//获取第i个备选服务器的流
			RequestStream<Request> const* thisStream = &alternatives->get(i, channel);
			if (!IFailureMonitor::failureMonitor().getState(thisStream->getEndpoint()).failed) {
				auto& qd = model->getMeasurement(thisStream->getEndpoint().token.first());
				if (now() > qd.failedUntil) {
					double thisMetric = qd.smoothOutstanding.smoothTotal();
					double thisTime = qd.latency;
					if (FLOW_KNOBS->LOAD_BALANCE_PENALTY_IS_BAD && qd.penalty > 1.001) {
						// When a server wants to penalize itself (the default
						// penalty value is 1.0), consider this server as bad.
						// penalty is sent from server.
						++badServers;
					}

					if (thisMetric < bestMetric) {
						if (i != bestAlt) {
							nextAlt = bestAlt;
							nextMetric = bestMetric;
							nextTime = bestTime;
						}
						bestAlt = i;
						bestMetric = thisMetric;
						bestTime = thisTime;
					} else if (thisMetric < nextMetric) {
						nextAlt = i;
						nextMetric = thisMetric;
						nextTime = thisTime;
					}
				} else {
					++badServers;
				}
			} else {
				++badServers;
			}
		}
		if (nextMetric > 1e8) {
			// If we still don't have a second best choice to issue request to,
			// go through all the remote servers again, since we may have
			// skipped it.
			for (int i = alternatives->countBest(); i < alternatives->size(); i++) {
				RequestStream<Request> const* thisStream = &alternatives->get(i, channel);
				if (!IFailureMonitor::failureMonitor().getState(thisStream->getEndpoint()).failed) {
					auto& qd = model->getMeasurement(thisStream->getEndpoint().token.first());
					if (now() > qd.failedUntil) {
						double thisMetric = qd.smoothOutstanding.smoothTotal();
						double thisTime = qd.latency;

						if (thisMetric < nextMetric) {
							nextAlt = i;
							nextMetric = thisMetric;
							nextTime = thisTime;
						}
					}
				}
			}
		}

		if (nextTime < 1e9) {
			// Decide when to send the request to the second best choice.
			if (bestTime > FLOW_KNOBS->INSTANT_SECOND_REQUEST_MULTIPLIER *
			                   (model->secondMultiplier * (nextTime) + FLOW_KNOBS->BASE_SECOND_REQUEST_TIME)) {
				secondDelay = Void();
			} else {
				secondDelay = delay(model->secondMultiplier * nextTime + FLOW_KNOBS->BASE_SECOND_REQUEST_TIME);
			}
		} else {
			secondDelay = Never();
		}
	}

	state int startAlt = nextAlt;
	state int startDistance = (bestAlt + alternatives->size() - startAlt) % alternatives->size();

	state int numAttempts = 0;
	state double backoff = 0;
	state bool triedAllOptions = false;
	// Issue requests to selected servers.
	//开始循环尝试发送请求
	loop {
		//如果请求持续时间超过阈值，记录一个警告事件，并详细记录尝试次数、退避时间等信息
		if (now() - startTime > (g_network->isSimulated() ? 30.0 : 600.0)) {
			TraceEvent ev(g_network->isSimulated() ? SevWarn : SevWarnAlways, "LoadBalanceTooLong");
			ev.suppressFor(1.0);
			ev.detail("Duration", now() - startTime);
			ev.detail("NumAttempts", numAttempts);
			ev.detail("Backoff", backoff);
			ev.detail("TriedAllOptions", triedAllOptions);
			if (ev.isEnabled()) {
				ev.log();
				for (int alternativeNum = 0; alternativeNum < alternatives->size(); alternativeNum++) {
					RequestStream<Request> const* thisStream = &alternatives->get(alternativeNum, channel);
					TraceEvent(SevWarn, "LoadBalanceTooLongEndpoint")
					    .detail("Addr", thisStream->getEndpoint().getPrimaryAddress())
					    .detail("Token", thisStream->getEndpoint().token)
					    .detail("Failed", IFailureMonitor::failureMonitor().getState(thisStream->getEndpoint()).failed);
				}
			}
		}

		// Find an alternative, if any, that is not failed, starting with
		// nextAlt. This logic matters only if model == nullptr. Otherwise, the
		// bestAlt and nextAlt have been decided.
		//尝试找到一个可用的服务器
		state RequestStream<Request> const* stream = NULL;
		for (int alternativeNum = 0; alternativeNum < alternatives->size(); alternativeNum++) {
			int useAlt = nextAlt;
			if (nextAlt == startAlt)
				useAlt = bestAlt;
			else if ((nextAlt + alternatives->size() - startAlt) % alternatives->size() <= startDistance)
				useAlt = (nextAlt + alternatives->size() - 1) % alternatives->size();

			stream = &alternatives->get(useAlt, channel);
			if (!IFailureMonitor::failureMonitor().getState(stream->getEndpoint()).failed &&
			    (!firstRequestEndpoint.present() || stream->getEndpoint().token.first() != firstRequestEndpoint.get()))
				break;
			nextAlt = (nextAlt + 1) % alternatives->size();
			if (nextAlt == startAlt)
				triedAllOptions = true;
			stream = NULL;
		}
		//如果所有服务器都失败了，等待其中一个恢复，根据不同的延迟策略，决定等待多长时间或立刻抛出异常
		if (!stream && !firstRequestData.isValid()) {
			// Everything is down!  Wait for someone to be up.

			vector<Future<Void>> ok(alternatives->size());
			for (int i = 0; i < ok.size(); i++) {
				ok[i] = IFailureMonitor::failureMonitor().onStateEqual(alternatives->get(i, channel).getEndpoint(),
				                                                       FailureStatus(false));
			}

			if (!alternatives->alwaysFresh()) {
				if (now() - g_network->networkInfo.newestAlternativesFailure >
				    FLOW_KNOBS->ALTERNATIVES_FAILURE_RESET_TIME) {
					g_network->networkInfo.oldestAlternativesFailure = now();
				}

				double delay = FLOW_KNOBS->ALTERNATIVES_FAILURE_MIN_DELAY;
				if (now() - g_network->networkInfo.lastAlternativesFailureSkipDelay >
				    FLOW_KNOBS->ALTERNATIVES_FAILURE_SKIP_DELAY) {
					g_network->networkInfo.lastAlternativesFailureSkipDelay = now();
				} else {
					double elapsed = now() - g_network->networkInfo.oldestAlternativesFailure;
					delay = std::max(delay,
					                 std::min(elapsed * FLOW_KNOBS->ALTERNATIVES_FAILURE_DELAY_RATIO,
					                          FLOW_KNOBS->ALTERNATIVES_FAILURE_MAX_DELAY));
					delay = std::max(delay,
					                 std::min(elapsed * FLOW_KNOBS->ALTERNATIVES_FAILURE_SLOW_DELAY_RATIO,
					                          FLOW_KNOBS->ALTERNATIVES_FAILURE_SLOW_MAX_DELAY));
				}

				// Making this SevWarn means a lot of clutter
				if (now() - g_network->networkInfo.newestAlternativesFailure > 1 ||
				    deterministicRandom()->random01() < 0.01) {
					TraceEvent("AllAlternativesFailed")
					    .detail("Interval", FLOW_KNOBS->CACHE_REFRESH_INTERVAL_WHEN_ALL_ALTERNATIVES_FAILED)
					    .detail("Alternatives", alternatives->description())
					    .detail("Delay", delay);
				}

				g_network->networkInfo.newestAlternativesFailure = now();

				choose {
					when(wait(quorum(ok, 1))) {}
					when(wait(::delayJittered(delay))) { throw all_alternatives_failed(); }
				}
			} else {
				wait(quorum(ok, 1));
			}

			numAttempts = 0; // now that we've got a server back, reset the backoff
		//如果只有第一个请求可用，处理第一个请求的响应
		} else if (!stream) {
			// Only the first location is available.
			ErrorOr<REPLY_TYPE(Request)> result = wait(firstRequestData.response);
			if (firstRequestData.checkAndProcessResult(atMostOnce)) {
				return result.get();
			}

			firstRequestEndpoint = Optional<uint64_t>();
		//如果第一个请求已发送且响应时间过长，发送第二个请求，等待其中之一响应
		} else if (firstRequestData.isValid()) {
			// Issue a second request, the first one is taking a long time.
			secondRequestData.startRequest(backoff, triedAllOptions, stream, request, model);
			state bool firstFinished = false;

			loop choose {
				when(ErrorOr<REPLY_TYPE(Request)> result =
				         wait(firstRequestData.response.isValid() ? firstRequestData.response : Never())) {
					if (firstRequestData.checkAndProcessResult(atMostOnce)) {
						return result.get();
					}

					firstRequestEndpoint = Optional<uint64_t>();
					firstFinished = true;
				}
				//等待请求的响应response
				when(ErrorOr<REPLY_TYPE(Request)> result = wait(secondRequestData.response)) {
					if (secondRequestData.checkAndProcessResult(atMostOnce)) {
						return result.get();
					}

					break;
				}
			}

			if (++numAttempts >= alternatives->size()) {
				backoff = std::min(
				    FLOW_KNOBS->LOAD_BALANCE_MAX_BACKOFF,
				    std::max(FLOW_KNOBS->LOAD_BALANCE_START_BACKOFF, backoff * FLOW_KNOBS->LOAD_BALANCE_BACKOFF_RATE));
			}
		} else {
			// Issue a request, if it takes too long to get a reply, go around the loop
			//这里是重点，发送请求
			firstRequestData.startRequest(backoff, triedAllOptions, stream, request, model);
			firstRequestEndpoint = stream->getEndpoint().token.first();

			loop {
				choose {
					when(ErrorOr<REPLY_TYPE(Request)> result = wait(firstRequestData.response)) {
						if (model) {
							model->secondMultiplier =
							    std::max(model->secondMultiplier - FLOW_KNOBS->SECOND_REQUEST_MULTIPLIER_DECAY, 1.0);
							model->secondBudget =
							    std::min(model->secondBudget + FLOW_KNOBS->SECOND_REQUEST_BUDGET_GROWTH,
							             FLOW_KNOBS->SECOND_REQUEST_MAX_BUDGET);
						}

						if (firstRequestData.checkAndProcessResult(atMostOnce)) {
							return result.get();
						}

						firstRequestEndpoint = Optional<uint64_t>();
						break;
					}
					when(wait(secondDelay)) {
						secondDelay = Never();
						if (model && model->secondBudget >= 1.0) {
							model->secondMultiplier += FLOW_KNOBS->SECOND_REQUEST_MULTIPLIER_GROWTH;
							model->secondBudget -= 1.0;
							break;
						}
					}
				}
			}

			if (++numAttempts >= alternatives->size()) {
				backoff = std::min(
				    FLOW_KNOBS->LOAD_BALANCE_MAX_BACKOFF,
				    std::max(FLOW_KNOBS->LOAD_BALANCE_START_BACKOFF, backoff * FLOW_KNOBS->LOAD_BALANCE_BACKOFF_RATE));
			}
		}
		//更新备选服务器索引并重置相关状态
		nextAlt = (nextAlt + 1) % alternatives->size();
		if (nextAlt == startAlt)
			triedAllOptions = true;
		resetReply(request, taskID);
		secondDelay = Never();
	}
}

// Subclasses must initialize all members in their default constructors
// Subclasses must serialize all members
struct BasicLoadBalancedReply {
	int processBusyTime;
	BasicLoadBalancedReply() : processBusyTime(0) {}
};

Optional<BasicLoadBalancedReply> getBasicLoadBalancedReply(const BasicLoadBalancedReply* reply);
Optional<BasicLoadBalancedReply> getBasicLoadBalancedReply(const void*);

// A simpler version of LoadBalance that does not send second requests where the list of servers are always fresh
//在多个可能的服务器之间进行负载均衡，以处理传入的请求。该函数没有发送第二个请求，并假设服务器列表总是最新的
//实现了一个基本的负载均衡算法，通过循环遍历多个服务器，找到一个可用的服务器来处理请求。如果所有服务器都不可用，等待至少一个服务器恢复可用。通过backoff机制和错误处理，提高了请求的成功率和稳定性
ACTOR template <class Interface, class Request, class Multi>
Future<REPLY_TYPE(Request)> basicLoadBalance(Reference<ModelInterface<Multi>> alternatives,//可用的服务器接口引用
                                             RequestStream<Request> Interface::*channel,//接口的请求流指针
                                             Request request = Request(),//要发送的请求
                                             TaskPriority taskID = TaskPriority::DefaultPromiseEndpoint,//任务优先级
                                             bool atMostOnce = false) {
	setReplyPriority(request, taskID);
	if (!alternatives)
		return Never();

	ASSERT(alternatives->size() && alternatives->alwaysFresh());

	state int bestAlt = alternatives->getBest(); //最佳服务器索引
	//下一个尝试的服务器索引
	state int nextAlt = deterministicRandom()->randomInt(0, std::max(alternatives->size() - 1, 1));
	if (nextAlt >= bestAlt)
		nextAlt++;
	//初始尝试的服务器索引
	state int startAlt = nextAlt;
	//从最佳服务器到初始服务器的距离
	state int startDistance = (bestAlt + alternatives->size() - startAlt) % alternatives->size();
	state int numAttempts = 0;//尝试次数
	state double backoff = 0;//回避时间
	state int useAlt;//当前使用的服务器索引
	loop {
		// Find an alternative, if any, that is not failed, starting with nextAlt
		//从nextAlt开始，尝试找到一个可用的服务器
		state RequestStream<Request> const* stream = NULL;
		for (int alternativeNum = 0; alternativeNum < alternatives->size(); alternativeNum++) {
			useAlt = nextAlt;
			if (nextAlt == startAlt)
				useAlt = bestAlt;
			else if ((nextAlt + alternatives->size() - startAlt) % alternatives->size() <= startDistance)
				useAlt = (nextAlt + alternatives->size() - 1) % alternatives->size();
			//找到可用的服务器后，将其请求流赋值给stream并break
			stream = &alternatives->get(useAlt, channel);
			if (!IFailureMonitor::failureMonitor().getState(stream->getEndpoint()).failed)
				break;
			nextAlt = (nextAlt + 1) % alternatives->size();
			stream = NULL;
		}
		//如果所有服务器都不可用，等待至少一个服务器恢复可用
		if (!stream) {
			// Everything is down!  Wait for someone to be up.

			vector<Future<Void>> ok(alternatives->size());
			for (int i = 0; i < ok.size(); i++) {
				ok[i] = IFailureMonitor::failureMonitor().onStateEqual(alternatives->get(i, channel).getEndpoint(),
				                                                       FailureStatus(false));
			}
			wait(quorum(ok, 1));

			numAttempts = 0; // now that we've got a server back, reset the backoff
		} else {
			//等待backoff的时间
			if (backoff > 0.0) {
				wait(delay(backoff));
			}
			//发送请求并等待响应
			ErrorOr<REPLY_TYPE(Request)> result = wait(stream->tryGetReply(request));
			//如果成功获取响应，更新服务器的负载信息并返回响应
			if (result.present()) {
				Optional<BasicLoadBalancedReply> loadBalancedReply = getBasicLoadBalancedReply(&result.get());
				if (loadBalancedReply.present()) {
					alternatives->updateRecent(useAlt, loadBalancedReply.get().processBusyTime);
				}

				return result.get();
			}
			//错误处理
			if (result.getError().code() != error_code_broken_promise &&
			    result.getError().code() != error_code_request_maybe_delivered) {
				throw result.getError();
			}

			if (atMostOnce) {
				throw request_maybe_delivered();
			}
			//设置backoff
			if (++numAttempts >= alternatives->size()) {
				backoff = std::min(
				    FLOW_KNOBS->LOAD_BALANCE_MAX_BACKOFF,
				    std::max(FLOW_KNOBS->LOAD_BALANCE_START_BACKOFF, backoff * FLOW_KNOBS->LOAD_BALANCE_BACKOFF_RATE));
			}
		}
		//更新nextAlt为下一个服务器索引，重置请求的回复状态
		nextAlt = (nextAlt + 1) % alternatives->size();
		resetReply(request, taskID);
	}
}

#include "flow/unactorcompiler.h"

#endif
