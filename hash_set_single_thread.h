/**
 * @file hash_set.h
 * @brief HashSet Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @modified by Huang Yukun, disable openmp parallel
 * @version 1.0.0
 * @date 2011-08-24
 */

#ifndef __SINGLETHREAD_CONTAINER_HASH_SET_H_

#define __SINGLETHREAD_CONTAINER_HASH_SET_H_

#include "hash.h"
#include "functional.h"
#include "hash_table_single_thread.h"

#include <functional>


/**
 * @brief It is a parallel hash set which has similar interface as stl set.
 * It is implemented based on parallel hash table (HashTable).
 *
 * @tparam Value
 * @tparam HashFunc
 */
template <typename Value, typename HashFunc = Hash<Value>,
         typename EqualKey = std::equal_to<Value> >
class HashSet
{
public:
    typedef HashTable<Value, Value, HashFunc, Identity<Value>, EqualKey> hash_table_type;
    typedef HashSet<Value, HashFunc, EqualKey> hash_set_type;

    typedef typename hash_table_type::key_type key_type;
    typedef typename hash_table_type::value_type value_type;
    typedef typename hash_table_type::size_type size_type;
    typedef typename hash_table_type::difference_type difference_type;

    typedef typename hash_table_type::reference reference; 
    typedef typename hash_table_type::const_reference const_reference;
    typedef typename hash_table_type::pointer pointer;
    typedef typename hash_table_type::const_pointer const_pointer;

    typedef typename hash_table_type::hash_func_type hash_func_type;
    typedef typename hash_table_type::get_key_func_type get_key_func_type;
    typedef typename hash_table_type::key_equal_func_type key_equal_func_type;

    typedef typename hash_table_type::iterator iterator;
    typedef typename hash_table_type::const_iterator const_iterator;

    explicit HashSet(const hash_func_type &hash = hash_func_type(),
            const key_equal_func_type key_equal = key_equal_func_type())
        : hash_table_(hash, Identity<value_type>(), key_equal)
    {}

    HashSet(const hash_set_type &hash_set)
        : hash_table_(hash_set.hash_table_)
    {}

    const hash_set_type &operator = (const hash_set_type &hash_set)
    { hash_table_ = hash_set.hash_table_; return *this; }

    iterator begin() { return hash_table_.begin(); }
    const_iterator begin() const { return hash_table_.begin(); }
    iterator end() { return hash_table_.end(); }
    const_iterator end() const { return hash_table_.end(); }

    std::pair<iterator, bool> insert(const value_type &value)
    { return hash_table_.insert_unique(value); }

    iterator find(const value_type &value)
    { return hash_table_.find(value); }

    const_iterator find(const value_type &value) const
    { return hash_table_.find(value); }

    size_type remove(const value_type &value)
    { return hash_table_.remove(value); }

    template <typename Predicator>
    size_type remove_if(Predicator &predicator)
    { return hash_table_.remove_if(predicator); }

    template <typename UnaryProc>
    UnaryProc &for_each(UnaryProc &op)
    { return hash_table_.for_each(op); }

    const hash_func_type &hash_func() const
    { return hash_table_.hash_func(); }
    const key_equal_func_type &key_equal_func() const
    { return hash_table_.key_equal_func(); }

    void reserve(size_type capacity)
    { hash_table_.reserve(capacity); }

    void swap(hash_set_type &hash_set)
    { if (this != &hash_set) hash_table_.swap(hash_set.hash_table_); }

    size_type size() const { return hash_table_.size(); }
    bool empty() const { return hash_table_.empty(); }

    void clear()
    { hash_table_.clear(); }

private:
    hash_table_type hash_table_;
};

namespace std
{
template <typename Value, typename HashFunc, typename EqualKey>
inline void swap(HashSet<Value, HashFunc, EqualKey> &x, 
        HashSet<Value, HashFunc, EqualKey> &y)
{ x.swap(y); }
}

#endif

