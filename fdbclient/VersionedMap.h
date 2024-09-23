/*
 * VersionedMap.h
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

#ifndef FDBCLIENT_VERSIONEDMAP_H
#define FDBCLIENT_VERSIONEDMAP_H
#pragma once

#include "flow/flow.h"
#include "flow/IndexedSet.h"
#include "fdbclient/FDBTypes.h"
#include "flow/IRandom.h"
#include "fdbclient/VersionedMap.actor.h"

// PTree is a persistent balanced binary tree implementation. It is based on a treap as a way to guarantee O(1) space
// for node insertion (rotating is asymptotically cheap), but the constant factors are very large.
//
// Each node has three pointers - the first two are its left and right children, respectively, and the third can be set
// to point to a newer version of the node. This third pointer allows us to maintain persistence without full path
// copying, and is employed to achieve O(1) space node insertion.
//
// PTree also supports efficient finger searches.
//PTree是一种持久化的平衡二叉树实现，基于treap数据结构，能够保证每次插入节点所需的空间是O1
//关键特性：
//1.持久化：所有对树的修改操作不会破坏原有版本的数据结构
//2.空间效率：使用版本指针避免每次插入时的全路径复制，提高空间使用效率
//3.良好的平衡性和搜索效率
namespace PTreeImpl {

#pragma warning(disable : 4800)
/**
 * PTree每次更新不会改变原有节点，而是创建一个新节点，通过第三个指针实现持久化更新，这种设计避免了完全路径复制，实现了O1空间的节点插入
 * PTree基于Treap树堆，结合了二叉查找树和堆的性质。每个节点有一个随机生成的优先级，插入和删除操作通过旋转维持树的平衡
 * PTree支持指针搜索，通过PTreeFinger类实现快速的前向和后向遍历
 */
template <class T>
struct PTree : public ReferenceCounted<PTree<T>>, FastAllocated<PTree<T>>, NonCopyable {
	uint32_t priority;
	//point[0]和poiont[1]分别指向左孩子和右孩子，point[2]指向新版本的节点
	Reference<PTree> pointer[3];
	Version lastUpdateVersion;
	bool updated; //是否更新过子树
	bool replacedPointer; //上次更新的是哪个子树 左/右
	T data; //存放具体的数据

	Reference<PTree> child(bool which, Version at) const {
		if (updated && lastUpdateVersion <= at && which == replacedPointer)
			return pointer[2];
		else
			return pointer[which];
	}
	Reference<PTree> left(Version at) const { return child(false, at); }
	Reference<PTree> right(Version at) const { return child(true, at); }

	PTree(const T& data, Version ver) : data(data), lastUpdateVersion(ver), updated(false) {
		//随机生成优先级
		priority = deterministicRandom()->randomUInt32();
	}
	PTree(uint32_t pri, T const& data, Reference<PTree> const& left, Reference<PTree> const& right, Version ver)
	  : priority(pri), data(data), lastUpdateVersion(ver), updated(false) {
		pointer[0] = left;
		pointer[1] = right;
	}

private:
	PTree(PTree const&);
};

//用于表示和操作树节点路径的辅助类，用于提高在树结构中进行查找和操作的效率
template <class T>
class PTreeFinger {
	using PTreeFingerEntry = PTree<T> const*;
	// This finger size supports trees with up to exp(96/4.3) ~= 4,964,514,749 entries.
	// see also: check().
	static constexpr size_t N = 96;
	//存储路径中的节点指针，支持树的最大深度为96
	PTreeFingerEntry entries_[N];
	//当前路径中节点指针的数量，记录entries数组中实际使用的元素个数
	size_t size_ = 0;
	//用于在范围查询操作中保存边界的大小
	size_t bound_sz_ = 0;

public:
	PTreeFinger() {}

	// Explicit copy constructors ensure we copy the live values in entries_.
	PTreeFinger(PTreeFinger const& f) { *this = f; }
	PTreeFinger(PTreeFinger&& f) { *this = f; }

	PTreeFinger& operator=(PTreeFinger const& f) {
		size_ = f.size_;
		bound_sz_ = f.bound_sz_;
		std::copy(f.entries_, f.entries_ + size_, entries_);
		return *this;
	}

	PTreeFinger& operator=(PTreeFinger&& f) {
		size_ = std::exchange(f.size_, 0);
		bound_sz_ = f.bound_sz_;
		std::copy(f.entries_, f.entries_ + size_, entries_);
		return *this;
	}

	size_t size() const { return size_; }
	PTree<T> const* back() const { return entries_[size_ - 1]; }
	void pop_back() { size_--; }
	void clear() { size_ = 0; }
	PTree<T> const* operator[](size_t i) const { return entries_[i]; }

	void resize(size_t sz) {
		size_ = sz;
		ASSERT(size_ < N);
	}

	void push_back(PTree<T> const* node) {
		entries_[size_++] = { node };
		ASSERT(size_ < N);
	}

	void push_for_bound(PTree<T> const* node, bool less) {
		push_back(node);
		bound_sz_ = less ? size_ : bound_sz_;
	}

	// remove the end of the finger so that the last entry is less than the probe
	void trim_to_bound() { size_ = bound_sz_; }
};

template <class T>
static Reference<PTree<T>> update(Reference<PTree<T>> const& node, //要更新的节点
                                  bool which, //表示要更新的是左子树还是右子树
                                  Reference<PTree<T>> const& ptr, //新的子树引用
                                  Version at) {  //更新的版本号
	//如果新子树节点与当前节点相同，则无需更新
	if (ptr.getPtr() == node->child(which, at).getPtr() /* && node->replacedVersion <= at*/) {
		return node;
	}
	//如果当前节点的版本与目标版本相同
	if (node->lastUpdateVersion == at) {
		//&& (!node->updated || node->replacedPointer==which)) {
		//如果已经更新过且更新的子树不同
		//简单来说就重新复制一下node，用ptr替换其左子树或右子树，同时原node的point2不再使用，清空
		if (node->updated && node->replacedPointer != which) {
			// We are going to have to copy this node, but its aux pointer will never be used again
			// and should drop its reference count
			Reference<PTree<T>> r;
			if (which) //更新右子树
				r = Reference<PTree<T>>(new PTree<T>(node->priority, node->data, node->child(0, at), ptr, at));
			else //更新左子树
				r = Reference<PTree<T>>(new PTree<T>(node->priority, node->data, ptr, node->child(1, at), at));
			node->pointer[2].clear();
			return r;
		} else { //直接更新节点的子树指针
			if (node->updated)
				node->pointer[2] = ptr;
			else
				node->pointer[which] = ptr;
			return node;
		}
	}
	//当前节点已经被更新过
	if (node->updated) {
		if (which)
			return Reference<PTree<T>>(new PTree<T>(node->priority, node->data, node->child(0, at), ptr, at));
		else
			return Reference<PTree<T>>(new PTree<T>(node->priority, node->data, ptr, node->child(1, at), at));
	} else {
		node->lastUpdateVersion = at;
		node->replacedPointer = which;
		node->pointer[2] = ptr;
		node->updated = true;
		return node;
	}
}

template <class T, class X>
bool contains(const Reference<PTree<T>>& p, Version at, const X& x) {
	if (!p)
		return false;
	int cmp = compare(x, p->data);
	bool less = cmp < 0;
	if (cmp == 0)
		return true;
	//从子树中查找
	return contains(p->child(!less, at), at, x);
}

// TODO: Remove the number of invocations of operator<, and replace with something closer to memcmp.
// and same for upper_bound.
//找到第一个不小于值x的节点，将路径记录在PTreeFinger中
template <class T, class X>
void lower_bound(const Reference<PTree<T>>& p, Version at, const X& x, PTreeFinger<T>& f) {
	if (!p) {
		f.trim_to_bound();
		return;
	}
	int cmp = compare(x, p->data);
	bool less = cmp < 0;
	f.push_for_bound(p.getPtr(), less);
	if (cmp == 0)
		return;
	lower_bound(p->child(!less, at), at, x, f);
}
//找到第一个不大于值x的节点，将路径记录在PTreeFinger中
template <class T, class X>
void upper_bound(const Reference<PTree<T>>& p, Version at, const X& x, PTreeFinger<T>& f) {
	if (!p) {
		f.trim_to_bound();
		return;
	}
	bool less = x < p->data;
	f.push_for_bound(p.getPtr(), less);
	upper_bound(p->child(!less, at), at, x, f);
}

template <class T, bool forward>
void move(Version at, PTreeFinger<T>& f) {
	ASSERT(f.size());
	const PTree<T>* n;
	//获取路径的最后一个节点
	n = f.back();
	if (n->child(forward, at)) { //如果在指定方向上有子节点
		n = n->child(forward, at).getPtr();
		do {
			f.push_back(n); //将访问的节点推入路径
			n = n->child(!forward, at).getPtr();
		} while (n); //继续沿反方向移动，直到没有子节点为止
	} else {
		do {
			n = f.back(); //获取路径中的最后一个节点
			f.pop_back(); //从路径中移除最后一个节点
		} while (f.size() && f.back()->child(forward, at).getPtr() == n); //回溯路径，直到找到一个在指定方向上有子节点的节点
	}
}

template <class T, bool forward>
int halfMove(Version at, PTreeFinger<T>& f) {
	// Post: f[:return_value] is the finger that would have been returned by move<forward>(at,f), and
	// f[:original_length_of_f] is unmodified
	ASSERT(f.size());
	const PTree<T>* n;
	n = f.back();
	if (n->child(forward, at)) {
		n = n->child(forward, at).getPtr();
		do {
			f.push_back(n);
			n = n->child(!forward, at).getPtr();
		} while (n);
		return f.size();
	} else {
		int s = f.size();
		do {
			n = f[s - 1];
			--s;
		} while (s && f[s - 1]->child(forward, at).getPtr() == n);
		return s;
	}
}

template <class T>
void next(Version at, PTreeFinger<T>& f) {
	move<T, true>(at, f);
}

template <class T>
void previous(Version at, PTreeFinger<T>& f) {
	move<T, false>(at, f);
}

template <class T>
int halfNext(Version at, PTreeFinger<T>& f) {
	return halfMove<T, true>(at, f);
}

template <class T>
int halfPrevious(Version at, PTreeFinger<T>& f) {
	return halfMove<T, false>(at, f);
}

template <class T>
T get(PTreeFinger<T>& f) {
	ASSERT(f.size());
	return f.back()->data;
}

// Modifies p to point to a PTree with x inserted
template <class T>
void insert(Reference<PTree<T>>& p, Version at, const T& x) {
	if (!p) { //直接创建一个根节点
		p = Reference<PTree<T>>(new PTree<T>(x, at));
	} else {
		//确定插入的方向
		bool direction = !(x < p->data);
		//获取对应方向的子节点
		Reference<PTree<T>> child = p->child(direction, at);
		//递归插入子节点
		insert(child, at, x);
		//更新当前节点p的子节点指针
		p = update(p, direction, child, at);
		//如果新插入的子节点优先级高于当前节点，进行旋转操作保持平衡
		if (p->child(direction, at)->priority > p->priority)
			rotate(p, at, !direction);
	}
}

template <class T>
Reference<PTree<T>> firstNode(const Reference<PTree<T>>& p, Version at) {
	if (!p)
		ASSERT(false);
	if (!p->left(at))
		return p;
	return firstNode(p->left(at), at);
}

template <class T>
Reference<PTree<T>> lastNode(const Reference<PTree<T>>& p, Version at) {
	if (!p)
		ASSERT(false);
	if (!p->right(at))
		return p;
	return lastNode(p->right(at), at);
}

template <class T, bool last>
void firstOrLastFinger(const Reference<PTree<T>>& p, Version at, PTreeFinger<T>& f) {
	if (!p)
		return;
	f.push_back(p.getPtr());
	firstOrLastFinger<T, last>(p->child(last, at), at, f);
}

template <class T>
void first(const Reference<PTree<T>>& p, Version at, PTreeFinger<T>& f) {
	return firstOrLastFinger<T, false>(p, at, f);
}

template <class T>
void last(const Reference<PTree<T>>& p, Version at, PTreeFinger<T>& f) {
	return firstOrLastFinger<T, true>(p, at, f);
}

// modifies p to point to a PTree with the root of p removed
template <class T>
void removeRoot(Reference<PTree<T>>& p, Version at) {
	if (!p->right(at))
		p = p->left(at);
	else if (!p->left(at))
		p = p->right(at);
	else {
		bool direction = p->right(at)->priority < p->left(at)->priority;
		rotate(p, at, direction);
		Reference<PTree<T>> child = p->child(direction, at);
		removeRoot(child, at);
		p = update(p, direction, child, at);
	}
}

// changes p to point to a PTree with x removed
template <class T, class X>
void remove(Reference<PTree<T>>& p, Version at, const X& x) {
	if (!p)
		ASSERT(false); // attempt to remove item not present in PTree
	int cmp = compare(x, p->data);
	if (cmp < 0) {
		Reference<PTree<T>> child = p->child(0, at);
		remove(child, at, x);
		p = update(p, 0, child, at);
	} else if (cmp > 0) {
		Reference<PTree<T>> child = p->child(1, at);
		remove(child, at, x);
		p = update(p, 1, child, at);
	} else {
		removeRoot(p, at);
	}
}

template <class T, class X>
void remove(Reference<PTree<T>>& p, Version at, const X& begin, const X& end) {
	if (!p)
		return;
	int beginDir, endDir;
	int beginCmp = compare(begin, p->data);
	if (beginCmp < 0)
		beginDir = -1;
	else if (beginCmp > 0)
		beginDir = +1;
	else
		beginDir = 0;
	if (!(p->data < end))
		endDir = -1;
	else
		endDir = +1;

	if (beginDir == endDir) {
		Reference<PTree<T>> child = p->child(beginDir == +1, at);
		remove(child, at, begin, end);
		p = update(p, beginDir == +1, child, at);
	} else {
		if (beginDir == -1) {
			Reference<PTree<T>> left = p->child(0, at);
			removeBeyond(left, at, begin, 1);
			p = update(p, 0, left, at);
		}
		if (endDir == +1) {
			Reference<PTree<T>> right = p->child(1, at);
			removeBeyond(right, at, end, 0);
			p = update(p, 1, right, at);
		}
		if (beginDir < endDir)
			removeRoot(p, at);
	}
}

template <class T, class X>
void removeBeyond(Reference<PTree<T>>& p, Version at, const X& pivot, bool dir) {
	if (!p)
		return;

	if ((p->data < pivot) ^ dir) {
		p = p->child(!dir, at);
		removeBeyond(p, at, pivot, dir);
	} else {
		Reference<PTree<T>> child = p->child(dir, at);
		removeBeyond(child, at, pivot, dir);
		p = update(p, dir, child, at);
	}
}

/*template<class T, class X>
void remove(Reference<PTree<T>>& p, Version at, const X& begin, const X& end) {
    Reference<PTree<T>> left, center, right;
    split(p, begin, left, center, at);
    split(center, end, center, right, at);
    p = append(left, right, at);
}*/

// inputs a PTree with the root node potentially violating the heap property
// modifies p to point to a valid PTree
template <class T>
void demoteRoot(Reference<PTree<T>>& p, Version at) {
	if (!p)
		ASSERT(false);

	uint32_t priority[2];
	for (int i = 0; i < 2; i++)
		if (p->child(i, at))
			priority[i] = p->child(i, at)->priority;
		else
			priority[i] = 0;

	bool higherDirection = priority[1] > priority[0];

	if (priority[higherDirection] < p->priority)
		return;

	// else, child(higherDirection) is a greater priority than us and the other child...
	rotate(p, at, !higherDirection);
	Reference<PTree<T>> child = p->child(!higherDirection, at);
	demoteRoot(child, at);
	p = update(p, !higherDirection, child, at);
}

template <class T>
Reference<PTree<T>> append(const Reference<PTree<T>>& left, const Reference<PTree<T>>& right, Version at) {
	if (!left)
		return right;
	if (!right)
		return left;

	Reference<PTree<T>> r = Reference<PTree<T>>(new PTree<T>(lastNode(left, at)->data, at));
	if (EXPENSIVE_VALIDATION) {
		ASSERT(r->data < firstNode(right, at)->data);
	}
	Reference<PTree<T>> a = left;
	remove(a, at, r->data);

	r->pointer[0] = a;
	r->pointer[1] = right;
	demoteRoot(r, at);
	return r;
}

template <class T, class X>
void split(Reference<PTree<T>> p, const X& x, Reference<PTree<T>>& left, Reference<PTree<T>>& right, Version at) {
	if (!p) {
		left = Reference<PTree<T>>();
		right = Reference<PTree<T>>();
		return;
	}

	if (p->data < x) {
		left = p;
		Reference<PTree<T>> lr = left->right(at);
		split(lr, x, lr, right, at);
		left = update(left, 1, lr, at);
	} else {
		right = p;
		Reference<PTree<T>> rl = right->left(at);
		split(rl, x, left, rl, at);
		right = update(right, 0, rl, at);
	}
}

template <class T>
void rotate(Reference<PTree<T>>& p, Version at, bool right) {
	auto r = p->child(!right, at);

	auto n1 = r->child(!right, at);
	auto n2 = r->child(right, at);
	auto n3 = p->child(right, at);

	auto newC = update(p, !right, n2, at);
	newC = update(newC, right, n3, at);
	p = update(r, !right, n1, at);
	p = update(p, right, newC, at);
}

template <class T>
void printTree(const Reference<PTree<T>>& p, Version at, int depth = 0) {
	if (p->left(at))
		printTree(p->left(at), at, depth + 1);
	for (int i = 0; i < depth; i++)
		printf("  ");
	// printf(":%s\n", describe(p->data.value.first).c_str());
	printf(":%s\n", describe(p->data.key).c_str());
	if (p->right(at))
		printTree(p->right(at), at, depth + 1);
}

template <class T>
void printTreeDetails(const Reference<PTree<T>>& p, int depth = 0) {
	// printf("Node %p (depth %d): %s\n", p.getPtr(), depth, describe(p->data.value.first).c_str());
	printf("Node %p (depth %d): %s\n", p.getPtr(), depth, describe(p->data.key).c_str());
	printf("  Left: %p\n", p->pointer[0].getPtr());
	printf("  Right: %p\n", p->pointer[1].getPtr());
	// if (p->pointer[2])
	if (p->updated)
		printf("  Version %lld %s: %p\n",
		       p->lastUpdateVersion,
		       p->replacedPointer ? "Right" : "Left",
		       p->pointer[2].getPtr());
	for (int i = 0; i < 3; i++)
		if (p->pointer[i])
			printTreeDetails(p->pointer[i], depth + 1);
}

/*static int depth(const Reference<PTree<int>>& p, Version at) {
    if (!p) return 0;
    int d1 = depth(p->left(at), at) + 1;
    int d2 = depth(p->right(at), at) + 1;
    return d1 > d2 ? d1 : d2;
}*/

template <class T>
void validate(const Reference<PTree<T>>& p, Version at, T* min, T* max, int& count, int& height, int depth = 0) {
	if (!p) {
		height = 0;
		return;
	}
	ASSERT((!min || *min <= p->data) && (!max || p->data <= *max));
	for (int i = 0; i < 2; i++) {
		if (p->child(i, at))
			ASSERT(p->child(i, at)->priority <= p->priority);
	}

	++count;
	int h1, h2;
	validate(p->left(at), at, min, &p->data, count, h1, depth + 1);
	validate(p->right(at), at, &p->data, max, count, h2, depth + 1);
	height = std::max(h1, h2) + 1;
}

template <class T>
void check(const Reference<PTree<T>>& p) {
	int count = 0, height;
	validate(p, (T*)0, (T*)0, count, height);
	if (count && height > 4.3 * log(double(count))) {
		// printf("height %d; count %d\n", height, count);
		ASSERT(false);
	}
}

// Remove pointers to any child nodes that have been updated at or before the given version
// This essentially gets rid of node versions that will never be read (beyond 5s worth of versions)
// TODO look into making this per-version compaction. (We could keep track of updated nodes at each version for example)
template <class T>
void compact(Reference<PTree<T>>& p, Version newOldestVersion) {
	if (!p) {
		return;
	}
	if (p->updated && p->lastUpdateVersion <= newOldestVersion) {
		/* If the node has been updated, figure out which pointer was repalced. And delete that pointer */
		auto which = p->replacedPointer;
		p->pointer[which] = p->pointer[2];
		p->updated = false;
		p->pointer[2] = Reference<PTree<T>>();
		// p->pointer[which] = Reference<PTree<T>>();
	}
	Reference<PTree<T>> left = p->left(newOldestVersion);
	Reference<PTree<T>> right = p->right(newOldestVersion);
	compact(left, newOldestVersion);
	compact(right, newOldestVersion);
}

} // namespace PTreeImpl

class ValueOrClearToRef {
public:
	static ValueOrClearToRef value(ValueRef const& v) { return ValueOrClearToRef(v, false); }
	static ValueOrClearToRef clearTo(KeyRef const& k) { return ValueOrClearToRef(k, true); }

	bool isValue() const { return !isClear; };
	bool isClearTo() const { return isClear; }

	ValueRef const& getValue() const {
		ASSERT(isValue());
		return item;
	};
	KeyRef const& getEndKey() const {
		ASSERT(isClearTo());
		return item;
	};

private:
	ValueOrClearToRef(StringRef item, bool isClear) : item(item), isClear(isClear) {}
	StringRef item;
	bool isClear;
};



/**
 * VersionedMap采用PTree来存储数据，每个版本对应一个root，通过维护不同版本的root，以实现版本化数据管理
 * PTree允许在插入和删除操作时不复制整个路径，从而在O1空间复杂度下实现节点插入
 * 
 */

// VersionedMap provides an interface to a partially persistent tree, allowing you to read the values at a particular
// version, create new versions, modify the current version of the tree, and forget versions prior to a specific
// version.
template <class K, class T>
class VersionedMap : NonCopyable {
	// private:
public:
	//持久平衡二叉树，用于存储键值对及其版本信息
	typedef PTreeImpl::PTree<MapPair<K, std::pair<T, Version>>> PTreeT;
	//用于高效遍历树的指针
	typedef PTreeImpl::PTreeFinger<MapPair<K, std::pair<T, Version>>> PTreeFingerT;
	typedef Reference<PTreeT> Tree;

	Version oldestVersion, latestVersion;

	// This deque keeps track of PTree root nodes at various versions. Since the
	// versions increase monotonically, the deque is implicitly sorted and hence
	// binary-searchable.
	//存储每个版本对应的PTree的根节点
	std::deque<std::pair<Version, Tree>> roots;

	struct rootsComparator {
		bool operator()(const std::pair<Version, Tree>& value, const Version& key) { return (value.first < key); }
		bool operator()(const Version& key, const std::pair<Version, Tree>& value) { return (key < value.first); }
	};

	//获取指定版本的PTree
	Tree const& getRoot(Version v) const {
		auto r = upper_bound(roots.begin(), roots.end(), v, rootsComparator());
		--r;
		return r->second;
	}

	// For each item in the versioned map, 4 PTree nodes are potentially allocated:
	static const int overheadPerItem = nextFastAllocatedSize(sizeof(PTreeT)) * 4;
	struct iterator;

	VersionedMap() : oldestVersion(0), latestVersion(0) { roots.emplace_back(0, Tree()); }
	VersionedMap(VersionedMap&& v) BOOST_NOEXCEPT : oldestVersion(v.oldestVersion),
	                                                latestVersion(v.latestVersion),
	                                                roots(std::move(v.roots)) {}
	void operator=(VersionedMap&& v) BOOST_NOEXCEPT {
		oldestVersion = v.oldestVersion;
		latestVersion = v.latestVersion;
		roots = std::move(v.roots);
	}

	Version getLatestVersion() const { return latestVersion; }
	Version getOldestVersion() const { return oldestVersion; }

	// front element should be the oldest version in the deque, hence the next oldest should be at index 1
	Version getNextOldestVersion() const { return roots[1]->first; }

	//删除早于指定版本的所有版本数据
	void forgetVersionsBefore(Version newOldestVersion) {
		ASSERT(newOldestVersion <= latestVersion);
		auto r = upper_bound(roots.begin(), roots.end(), newOldestVersion, rootsComparator());
		auto upper = r;
		--r;
		// if the specified newOldestVersion does not exist, insert a new
		// entry-pair with newOldestVersion and the root from next lower version
		if (r->first != newOldestVersion) {
			r = roots.emplace(upper, newOldestVersion, getRoot(newOldestVersion));
		}

		UNSTOPPABLE_ASSERT(r->first == newOldestVersion);
		roots.erase(roots.begin(), r);
		oldestVersion = newOldestVersion;
	}

	Future<Void> forgetVersionsBeforeAsync(Version newOldestVersion, TaskPriority taskID = TaskPriority::DefaultYield) {
		ASSERT(newOldestVersion <= latestVersion);
		auto r = upper_bound(roots.begin(), roots.end(), newOldestVersion, rootsComparator());
		auto upper = r;
		--r;
		// if the specified newOldestVersion does not exist, insert a new
		// entry-pair with newOldestVersion and the root from next lower version
		if (r->first != newOldestVersion) {
			r = roots.emplace(upper, newOldestVersion, getRoot(newOldestVersion));
		}

		UNSTOPPABLE_ASSERT(r->first == newOldestVersion);

		std::vector<Tree> toFree;
		toFree.reserve(10000);
		auto newBegin = r;
		Tree* lastRoot = nullptr;
		for (auto root = roots.begin(); root != newBegin; ++root) {
			if (root->second) {
				if (lastRoot != nullptr && root->second == *lastRoot) {
					(*lastRoot).clear();
				}
				if (root->second->isSoleOwner()) {
					toFree.push_back(root->second);
				}
				lastRoot = &root->second;
			}
		}

		roots.erase(roots.begin(), newBegin);
		oldestVersion = newOldestVersion;
		return deferredCleanupActor(toFree, taskID);
	}

public:
	//创建新的版本
	void createNewVersion(Version version) { // following sets and erases are into the given version, which may now be
		                                     // passed to at().  Must be called in monotonically increasing order.
		if (version > latestVersion) {
			latestVersion = version;
			Tree r = getRoot(version);
			roots.emplace_back(version, r);
		} else
			ASSERT(version == latestVersion);
	}

	// insert() and erase() invalidate atLatest() and all iterators into it
	void insert(const K& k, const T& t) { insert(k, t, latestVersion); }
	void insert(const K& k, const T& t, Version insertAt) {
		if (PTreeImpl::contains(roots.back().second, latestVersion, k))
			PTreeImpl::remove(roots.back().second,
			                  latestVersion,
			                  k); // FIXME: Make PTreeImpl::insert do this automatically  (see also WriteMap.h FIXME)
		PTreeImpl::insert(
		    roots.back().second, latestVersion, MapPair<K, std::pair<T, Version>>(k, std::make_pair(t, insertAt)));
	}
	void erase(const K& begin, const K& end) { PTreeImpl::remove(roots.back().second, latestVersion, begin, end); }
	void erase(const K& key) { // key must be present
		PTreeImpl::remove(roots.back().second, latestVersion, key);
	}
	void erase(iterator const& item) { // iterator must be in latest version!
		// SOMEDAY: Optimize to use item.finger and avoid repeated search
		K key = item.key();
		erase(key);
	}

	void printDetail() { PTreeImpl::printTreeDetails(roots.back().second, 0); }

	void printTree(Version at) { PTreeImpl::printTree(roots.back().second, at, 0); }

	//压缩旧版本数据，移除不必要的节点
	void compact(Version newOldestVersion) {
		ASSERT(newOldestVersion <= latestVersion);
		// auto newBegin = roots.lower_bound(newOldestVersion);
		auto newBegin = lower_bound(roots.begin(), roots.end(), newOldestVersion, rootsComparator());
		for (auto root = roots.begin(); root != newBegin; ++root) {
			if (root->second)
				PTreeImpl::compact(root->second, newOldestVersion);
		}
		// printf("\nPrinting the tree at latest version after compaction.\n");
		// PTreeImpl::printTreeDetails(roots.back().second(), 0);
	}

	// for(auto i = vm.at(version).lower_bound(range.begin); i < range.end; ++i)
	struct iterator {
		explicit iterator(Tree const& root, Version at) : root(root), at(at) {}

		K const& key() const { return finger.back()->data.key; }
		Version insertVersion() const {
			return finger.back()->data.value.second;
		} // Returns the version at which the current item was inserted
		operator bool() const { return finger.size() != 0; }
		bool operator<(const K& key) const { return this->key() < key; }

		T const& operator*() { return finger.back()->data.value.first; }
		T const* operator->() { return &finger.back()->data.value.first; }
		void operator++() {
			if (finger.size())
				PTreeImpl::next(at, finger);
			else
				PTreeImpl::first(root, at, finger);
		}
		void operator--() {
			if (finger.size())
				PTreeImpl::previous(at, finger);
			else
				PTreeImpl::last(root, at, finger);
		}
		bool operator==(const iterator& r) const {
			if (finger.size() && r.finger.size())
				return finger.back() == r.finger.back();
			else
				return finger.size() == r.finger.size();
		}
		bool operator!=(const iterator& r) const {
			if (finger.size() && r.finger.size())
				return finger.back() != r.finger.back();
			else
				return finger.size() != r.finger.size();
		}

	private:
		friend class VersionedMap<K, T>;
		Tree root;
		Version at;
		PTreeFingerT finger;
	};

	class ViewAtVersion {
	public:
		ViewAtVersion(Tree const& root, Version at) : root(root), at(at) {}

		iterator begin() const {
			iterator i(root, at);
			PTreeImpl::first(root, at, i.finger);
			return i;
		}
		iterator end() const { return iterator(root, at); }

		// Returns x such that key==*x, or end()
		template <class X>
		iterator find(const X& key) const {
			iterator i(root, at);
			PTreeImpl::lower_bound(root, at, key, i.finger);
			if (i && i.key() == key)
				return i;
			else
				return end();
		}

		// Returns the smallest x such that *x>=key, or end()
		template <class X>
		iterator lower_bound(const X& key) const {
			iterator i(root, at);
			PTreeImpl::lower_bound(root, at, key, i.finger);
			return i;
		}

		// Returns the smallest x such that *x>key, or end()
		template <class X>
		iterator upper_bound(const X& key) const {
			iterator i(root, at);
			PTreeImpl::upper_bound(root, at, key, i.finger);
			return i;
		}

		// Returns the largest x such that *x<=key, or end()
		template <class X>
		iterator lastLessOrEqual(const X& key) const {
			iterator i(root, at);
			PTreeImpl::upper_bound(root, at, key, i.finger);
			--i;
			return i;
		}

		// Returns the largest x such that *x<key, or end()
		template <class X>
		iterator lastLess(const X& key) const {
			iterator i(root, at);
			PTreeImpl::lower_bound(root, at, key, i.finger);
			--i;
			return i;
		}

		void validate() {
			int count = 0, height = 0;
			PTreeImpl::validate<MapPair<K, std::pair<T, Version>>>(root, at, NULL, NULL, count, height);
			if (height > 100)
				TraceEvent(SevWarnAlways, "DiabolicalPTreeSize").detail("Size", count).detail("Height", height);
		}
		

	private:
		Tree root;
		Version at;
	};

	ViewAtVersion at(Version v) const { return ViewAtVersion(getRoot(v), v); }
	ViewAtVersion atLatest() const { return ViewAtVersion(roots.back().second, latestVersion); }

	bool isClearContaining(ViewAtVersion const& view, KeyRef key) {
		auto i = view.lastLessOrEqual(key);
		return i && i->isClearTo() && i->getEndKey() > key;
	}

	std::vector<std::pair<Version, T>> getAllVersionsOfKey(const K& key, Version version) {
		std::vector<std::pair<Version, T>> result;

		// 使用 lower_bound 查找大于等于指定版本的第一个版本根节点
		auto r = std::lower_bound(roots.begin(), roots.end(), version, rootsComparator());
		
		// 从找到的版本开始，遍历所有版本
		for (; r != roots.end(); ++r) {
			Version currentVersion = r->first;
			Tree currentRoot = r->second;

			// 使用 PTreeFinger 来查找指定版本中的键
			PTreeFingerT finger;
			PTreeImpl::lower_bound(currentRoot, currentVersion, key, finger);

			// 检查找到的键是否与目标键相等
			if (finger.size() && finger.back()->data.key == key) {
				// 记录当前版本和对应的值
				result.emplace_back(currentVersion, finger.back()->data.value.first);
			}
		}

		return result; // 返回所有找到的版本和值
	}

	// TODO: getHistory?
};

#endif
