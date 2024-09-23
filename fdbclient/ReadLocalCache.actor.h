#include <unordered_map>
#include <chrono>
#include <iostream>
#include <fdbclient/FDBTypes.h>
#include <fdbclient/ReadYourWrites.h>

// 使用 FoundationDB 的命名空间
using namespace fdb;

// 缓存条目结构体
struct CacheEntry {
    Value value;
    Version version;
    std::chrono::steady_clock::time_point timestamp;
};

// 本地缓存类，支持存储同一 key 的多个版本，带有 LRU 淘汰策略
class LocalCache {
public:
    LocalCache(size_t maxSize, std::chrono::seconds ttl)
        : maxSize(maxSize), ttl(ttl) {}

    // 获取特定 key 和版本的值
    Optional<Value> get(const Key& key, Version version) {
        auto it = cache.find(key);
        if (it != cache.end()) {
            auto& versionedCache = it->second;

            auto versionIt = versionedCache.find(version);
            if (versionIt != versionedCache.end()) {
                // 检查缓存条目是否过期
                if (std::chrono::steady_clock::now() - versionIt->second.timestamp < ttl) {
                    // 将此 key 移动到最近使用的位置
                    touch(key);
                    return versionIt->second.value;
                } else {
                    versionedCache.erase(versionIt);  // 移除过期条目
                }
            }
        }
        return {};
    }

    // 存储特定 key 和版本的值
    void put(const Key& key, const Value& value, Version version) {
        // 如果缓存大小超限，进行淘汰
        if (cache.size() >= maxSize) {
            evict();
        }

        // 存储 key 对应的版本数据，并更新 LRU 顺序
        cache[key][version] = { value, version, std::chrono::steady_clock::now() };
        touch(key);
    }

private:
    // 访问或插入 key 时，更新其在 LRU 顺序中的位置
    void touch(const Key& key) {
        // 如果 key 已存在于 LRU 列表中，删除旧的位置
        if (keyPositions.find(key) != keyPositions.end()) {
            lruList.erase(keyPositions[key]);
        }
        // 将 key 插入到 LRU 列表的头部
        lruList.push_front(key);
        keyPositions[key] = lruList.begin();  // 更新 key 在 LRU 列表中的位置
    }

    // LRU 淘汰策略：移除最久未使用的 key
    void evict() {
        if (!lruList.empty()) {
            // 从 LRU 列表的尾部获取最久未使用的 key
            Key keyToEvict = lruList.back();
            lruList.pop_back();  // 移除 LRU 列表中的该 key
            cache.erase(keyToEvict);  // 删除该 key 的所有版本数据
            keyPositions.erase(keyToEvict);  // 移除 key 位置的记录
        }
    }

    size_t maxSize;  // 最大缓存大小（按 key 计数）
    std::chrono::seconds ttl;  // 缓存时间
    std::unordered_map<Key, std::unordered_map<Version, CacheEntry>> cache;  // 多版本缓存
    std::list<Key> lruList;  // LRU 顺序记录
    std::unordered_map<Key, std::list<Key>::iterator> keyPositions;  // 记录每个 key 在 LRU 列表中的位置
};
