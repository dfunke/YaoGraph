#pragma once

#include <queue>
#include <utility>

template<class T, class Comp = std::less<T>>
class PriQueueD {
public:
    PriQueueD(size_t capacity, size_t log_degree = 3)
        : _capacity(capacity), _ldeg(log_degree) {
        _heap_element_vec.reserve(capacity);
        _heap_hindex_vec.reserve(capacity);
        _handle_vec.reserve(capacity);
        _free_handle_vec.reserve(capacity);
    }
    ~PriQueueD() {}

private:
    // identifier object created on push, and deleted on pop

public:
    using handle = size_t;

    const T &top() const {
        ASSERT(!empty());
        return _heap_element_vec[0];
    }

    bool empty() const { return _heap_element_vec.empty(); }
    size_t size() const { return _heap_element_vec.size(); }

    handle push(T value) {
        auto h = get_handle();

        auto pos = _heap_element_vec.size();
        // * First Iteration of sift_up ****************************************
        size_t par = parent(pos);
        if (_heap_element_vec.empty() || compare(_heap_element_vec[par], value)) {
            _heap_element_vec.emplace_back(std::move(value));
            _heap_hindex_vec.emplace_back(h);
            _handle_vec[h] = pos;
            return h;
        }
        _heap_element_vec.emplace_back(std::move(_heap_element_vec[par]));
        _heap_hindex_vec.emplace_back(_heap_hindex_vec[par]);
        _handle_vec[_heap_hindex_vec[par]] = pos;

        pos = sift_up(par, value);
        _heap_element_vec[pos] = std::move(value);
        _heap_hindex_vec[pos] = h;
        _handle_vec[_heap_hindex_vec[pos]] = pos;

        return _heap_hindex_vec[pos];
    }

    void pop() {
        ASSERT(!empty());

        auto temp_element = std::move(_heap_element_vec[_heap_element_vec.size() - 1]);
        auto temp_handle = _heap_hindex_vec[_heap_hindex_vec.size() - 1];

        _heap_element_vec.pop_back();
        _heap_hindex_vec.pop_back();

        _free_handle_vec.push_back(_heap_hindex_vec[0]);
        if (!_heap_element_vec.size()) return;

        auto pos = sift_down(0, temp_element);
        _heap_element_vec[pos] = std::move(temp_element);
        _heap_hindex_vec[pos] = temp_handle;
        _handle_vec[temp_handle] = pos;
    }

    void remove(handle h) {
        ASSERT(!empty());

        auto pos = _handle_vec[h];

        auto temp_element = std::move(_heap_element_vec[_heap_element_vec.size() - 1]);
        auto temp_handle = _heap_hindex_vec[_heap_hindex_vec.size() - 1];

        _heap_element_vec.pop_back();
        _heap_hindex_vec.pop_back();

        _free_handle_vec.push_back(h);
        if (pos >= _heap_element_vec.size()) return;

        if (compare(temp_element, _heap_element_vec[pos])) {
            pos = sift_up(pos, temp_element);
        } else {
            pos = sift_down(pos, temp_element);
        }
        _heap_element_vec[pos] = std::move(temp_element);
        _heap_hindex_vec[pos] = temp_handle;
        _handle_vec[temp_handle] = pos;
    }

    void change_key(handle h, T newvalue) {

        auto pos = _handle_vec[h];
        if (compare(newvalue, _heap_element_vec[pos])) {
            pos = sift_up(pos, newvalue);
        } else {
            pos = sift_down(pos, newvalue);
        }
        _heap_element_vec[pos] = std::move(newvalue);
        _heap_hindex_vec[pos] = h;
        _handle_vec[h] = pos;
    }

    const T &get_key(handle h) const {
        return _heap_element_vec[_handle_vec[h]];
    }

    bool verify() const {
        return verify_node(0);
    }

private:
    /* member definitions *****************************************************/
    //std::priority_queue<T, std::vector<T>, Comp> pq;
    size_t _capacity;
    size_t _ldeg;
    Comp _comp;
    std::vector<size_t> _heap_hindex_vec;
    std::vector<T> _heap_element_vec;
    std::vector<size_t> _handle_vec;
    std::vector<size_t> _free_handle_vec;


    bool compare(const T &v1, const T &v2) const {
        return _comp(v2, v1);
    }

    void move(size_t from, size_t too) {
        _heap_hindex_vec[too] = _heap_hindex_vec[from];
        _heap_element_vec[too] = std::move(_heap_element_vec[from]);
        _handle_vec[_heap_hindex_vec[too]] = too;
    }

    size_t get_handle() {
        if (_free_handle_vec.size()) {
            size_t b = _free_handle_vec[_free_handle_vec.size() - 1];
            _free_handle_vec.pop_back();
            return b;
        } else {
            if (_handle_vec.size() == _capacity) {
                _capacity *= 2;

                _heap_element_vec.reserve(_capacity);
                _heap_hindex_vec.reserve(_capacity);
                _handle_vec.reserve(_capacity);
            }

            _handle_vec.push_back(0ull);
            return _handle_vec.size() - 1;
        }
    }

    size_t parent(size_t index) const {
        return (index - 1) >> _ldeg;
    }

    size_t child(size_t index) const {
        return (index << _ldeg) + 1;
    }

    size_t sift_up(size_t pos, const T &val) {
        while (pos != 0) {
            size_t par = parent(pos);
            if (compare(_heap_element_vec[par], val)) return pos;
            move(par, pos);
            pos = par;
        }
        return 0;
    }

    size_t sift_down(size_t pos, const T &val) {
        // while not in leaf
        while (child(pos) < _heap_element_vec.size()) {
            size_t min_idx = child(pos);

            // find min child
            for (size_t i = child(pos) + 1;
                 i < child(pos + 1) && i < _heap_element_vec.size();
                 ++i) {
                min_idx = compare(_heap_element_vec[min_idx], _heap_element_vec[i]) ? min_idx : i;
            }

            if (compare(val, _heap_element_vec[min_idx])) return pos;
            //_vec[pos] = std::move(_vec[min_idx]);
            move(min_idx, pos);
            pos = min_idx;
        }
        return pos;
    }

    bool verify_node(size_t pos) const {
        // check all children
        for (size_t i = child(pos);
             i < child(pos + 1) && i < _heap_element_vec.size();
             ++i) {
            if (!compare(_heap_element_vec[pos], _heap_element_vec[i]) || !verify_node(i))
                return false;
        }

        return true;
    }
};


template<class T, class Comp = std::less<T>>
class PriQueueAdapter {

public:
    using pqType = PriQueueD<T, Comp>;
    using arrType = std::vector<T>;
    using handle = typename pqType::handle;

    PriQueueAdapter(arrType &&static_elements, size_t capacityPQ, size_t log_degree = 3)
        : arr(std::move(static_elements)), pq(capacityPQ, log_degree) {

        // initially sort the static_elements according to Comp
        std::sort(arr.begin(), arr.end(), comp);
        arrIdx = arr.size();
    }

public:
    const T &top() const {
        ASSERT(!empty());

        if (arrIdx > 0 && !pq.empty()) {
            // both are not empty, we need to compare

            if (compare(arr[arrIdx - 1], pq.top())) {
                return arr[arrIdx - 1];
            } else {
                return pq.top();
            }
        } else {
            // one is empty, can't be both per ASSERT

            if (!pq.empty()) {
                ASSERT(arrIdx == 0);
                return pq.top();
            } else {
                ASSERT(arrIdx > 0);
                return arr[arrIdx - 1];
            }
        }
    }

    bool empty() const { return arrIdx == 0 && pq.empty(); }
    size_t size() const { return arrIdx + pq.size(); }

    handle push(T value) {
        return pq.push(value);
    }

    void pop() {
        ASSERT(!empty());

        if (arrIdx > 0 && !pq.empty()) {
            // both are not empty, we need to compare

            if (compare(arr[arrIdx - 1], pq.top())) {
                arrIdx--;
            } else {
                pq.pop();
            }
        } else {
            // one is empty, can't be both per ASSERT

            if (!pq.empty()) {
                ASSERT(arrIdx == 0);
                pq.pop();
            } else {
                ASSERT(arrIdx > 0);
                arrIdx--;
            }
        }
    }

    void remove(handle h) {
        pq.remove(h);
    }

    void change_key(handle h, T newvalue) {
        pq.change_key(h, newvalue);
    }

    const T &get_key(handle h) const {
        return pq.get_key(h);
    }

    bool verify() const {
        return pq.verify();
    }

private:
    bool compare(const T &v1, const T &v2) const {
        return comp(v2, v1);
    }

private:
    arrType arr;
    pqType pq;
    typename arrType::size_type arrIdx;
    Comp comp;
};