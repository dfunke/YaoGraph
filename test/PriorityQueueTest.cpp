#include <gtest/gtest.h>

#include "utils/PriorityQueue.hpp"

#include "element_mock.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <string>

template<template<class, class> class Q>
struct D1 {
    template<class T, int kDegree, class Comp = std::less<>>
    static auto make(int size) {
        return Q<T, Comp>(size, kDegree);
    }
};

template<template<class, int, class> class Q>
struct D2 {
    template<class T, int kDegree, class Comp = std::less<>>
    static auto make(int size) {
        return Q<T, kDegree, Comp>(size);
    }
};

template<template<class, class> class Q>
constexpr D1<Q> D_Helper();

template<template<class, int, class> class Q>
constexpr D2<Q> D_Helper();

using D = decltype(D_Helper<PriQueueD>());

template<class P>
class PriorityQueueTest : public ::testing::Test {
protected:
    template<class T = int, int kDegree = 8, class Comp = std::less<>>
    auto make(int size = 300) const {
        return P::template make<T, kDegree, Comp>(size);
    }
};


using MyTypes = ::testing::Types<D>;
TYPED_TEST_CASE(PriorityQueueTest, MyTypes);

TYPED_TEST(PriorityQueueTest, ReportsSize) {
    auto q = this->make();
    ASSERT_TRUE(q.empty());
    ASSERT_EQ(q.size(), 0ul);

    q.push(0);
    ASSERT_FALSE(q.empty());
    ASSERT_EQ(q.size(), 1ul);

    q.push(1);
    q.pop();
    q.pop();
    ASSERT_TRUE(q.empty());
    ASSERT_EQ(q.size(), 0ul);

    for (unsigned int i = 1; i < 289; ++i) {
        q.push(i);
        ASSERT_EQ(q.size(), i);
    }
    ASSERT_FALSE(q.empty());
}

TYPED_TEST(PriorityQueueTest, ReturnsElementsInCorrectOrder) {
    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist;
    std::vector<int> values(300);
    auto q = this->make();

    for (int i = 0; i < 10; ++i) {
        std::generate(values.begin(), values.end(), [&] { return dist(rng); });

        for (auto &&v : values) q.push(v);
        ASSERT_EQ(q.size(), values.size());

        std::sort(values.begin(), values.end(), std::greater<>{});

        for (auto &&v : values) {
            ASSERT_EQ(q.top(), v);
            q.pop();
        }
        ASSERT_TRUE(q.empty());
    }
}

TYPED_TEST(PriorityQueueTest, CanInsertAfterPop) {
    const std::size_t n = 300;// should be even

    auto q = this->make();
    for (std::size_t i = 0; i < n; ++i)
        q.push(i);
    ASSERT_EQ(q.size(), n);

    for (auto i = n - 1; i >= n / 2; --i) {
        ASSERT_EQ(q.top(), static_cast<int>(i));
        q.pop();
    }
    ASSERT_EQ(q.size(), n / 2);

    for (auto i = 0ul; i < n / 2; ++i) {
        q.push(2 * n + i);
    }
    ASSERT_EQ(q.size(), n);

    for (auto i = n / 2; i > 0; --i) {
        ASSERT_EQ(q.top(), static_cast<int>(2 * n - 1 + i));
        q.pop();
    }

    for (auto i = n / 2; i > 0; --i) {
        ASSERT_EQ(q.top(), static_cast<int>(i - 1));
        q.pop();
    }
    ASSERT_TRUE(q.empty());
}

TYPED_TEST(PriorityQueueTest, WorksForNonTrivialDataTypes) {
    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist(0, 50);
    std::vector<std::string> values(128);
    auto q = this->template make<std::string>();

    for (int i = 0; i < 10; ++i) {
        std::generate(values.begin(), values.end(),
                      [&] { return std::string(dist(rng), 'x'); });

        for (auto &&v : values) q.push(v);
        ASSERT_EQ(q.size(), values.size());

        std::sort(values.begin(), values.end(), std::greater<>{});

        for (auto &&v : values) {
            ASSERT_EQ(q.top(), v);
            q.pop();
        }
        ASSERT_TRUE(q.empty());
    }
}

struct StringLengthComparator {
    bool operator()(const std::string &lhs, const std::string &rhs) const {
        if (lhs.size() == rhs.size())
            return lhs < rhs;
        return lhs.size() < rhs.size();
    }
    static bool reverse(const std::string &lhs, const std::string &rhs) {
        if (lhs.size() == rhs.size())
            return lhs > rhs;
        return lhs.size() > rhs.size();
    }
};

TYPED_TEST(PriorityQueueTest, WorksForDifferentComparator) {
    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist(0, 50);
    std::uniform_int_distribution<char> cdist('a', 'z');
    std::vector<std::string> values(128);
    auto q = this->template make<std::string, 8, StringLengthComparator>();

    for (int i = 0; i < 3; ++i) {
        std::generate(values.begin(), values.end(), [&] {
            std::string s(dist(rng), 'x');
            std::generate(s.begin(), s.end(), [&] { return cdist(rng); });
            return s;
        });

        for (auto &&v : values) q.push(v);
        ASSERT_EQ(q.size(), values.size());

        std::sort(values.begin(), values.end(), StringLengthComparator::reverse);

        for (auto &&v : values) {
            ASSERT_EQ(q.top(), v);
            q.pop();
        }
        ASSERT_TRUE(q.empty());
    }
}

TYPED_TEST(PriorityQueueTest, WorksForDifferentDegrees) {
    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist;
    std::vector<int> values(300);
    {
        auto q = this->template make<int, 3>();

        for (int i = 0; i < 5; ++i) {
            std::generate(values.begin(), values.end(), [&] { return dist(rng); });

            for (auto &&v : values) q.push(v);
            ASSERT_EQ(q.size(), values.size());

            std::sort(values.begin(), values.end(), std::greater<>{});

            for (auto &&v : values) {
                ASSERT_EQ(q.top(), v);
                q.pop();
            }
            ASSERT_TRUE(q.empty());
        }
    }
    {
        auto q = this->template make<int, 15>();

        for (int i = 0; i < 5; ++i) {
            std::generate(values.begin(), values.end(), [&] { return dist(rng); });

            for (auto &&v : values) q.push(v);
            ASSERT_EQ(q.size(), values.size());

            std::sort(values.begin(), values.end(), std::greater<>{});

            for (auto &&v : values) {
                ASSERT_EQ(q.top(), v);
                q.pop();
            }
            ASSERT_TRUE(q.empty());
        }
    }
}

TYPED_TEST(PriorityQueueTest, OPTIONAL_OnlyMoves) {
    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist;
    std::vector<ElementMock> values(300);
    std::vector<ElementMock> values_copy(300);
    auto q = this->template make<ElementMock>();

    for (int i = 0; i < 10; ++i) {
        std::generate(values.begin(), values.end(), [&] { return dist(rng); });
        values_copy = values;
        std::sort(values_copy.begin(), values_copy.end(), std::greater<>{});

        ElementMock::resetCounters();

        for (auto &&v : values) q.push(std::move(v));
        ASSERT_EQ(q.size(), values.size());
        ASSERT_EQ(ElementMock::copies, 0);
        EXPECT_EQ(ElementMock::default_constructions, 0);

        for (auto &&v : values_copy) {
            ASSERT_EQ(q.top(), v);
            q.pop();
        }
        EXPECT_EQ(ElementMock::copies, 0);
        EXPECT_EQ(ElementMock::default_constructions, 0);
        EXPECT_TRUE(q.empty());
    }
}

TEST(PriorityQueueTest, CanIncreaseKey) {
    std::vector<int> values(300);
    std::iota(values.begin(), values.end(), 0);

    auto q = D::make<int, 8>(300);
    std::vector<decltype(q)::handle> handles;
    for (auto v : values) handles.push_back(q.push(v));
    ASSERT_EQ(q.size(), values.size());

    std::mt19937_64 rng;
    std::uniform_int_distribution<int> dist(500, 100000);
    std::vector<int> new_values(values.size());
    std::generate(new_values.begin(), new_values.end(), [&] { return dist(rng); });

    for (auto i = 0ul; i < new_values.size(); ++i) {
        ASSERT_EQ(q.get_key(handles[i]), static_cast<int>(i));
        q.change_key(handles[i], new_values[i]);
    }
    ASSERT_EQ(q.size(), values.size());

    std::sort(new_values.begin(), new_values.end(), std::greater<>{});

    for (auto v : new_values) {
        ASSERT_EQ(q.top(), v);
        q.pop();
    }
    ASSERT_TRUE(q.empty());
}

TEST(PriorityQueueTest, CanRemove) {
    std::vector<int> values(300);
    std::iota(values.begin(), values.end(), 0);

    auto q = D::make<int, 8>(300);
    std::vector<decltype(q)::handle> handles;
    for (auto v : values) handles.push_back(q.push(v));
    ASSERT_EQ(q.size(), values.size());

    for (auto i = 0ul; i < handles.size(); i += 2) {
        ASSERT_EQ(q.get_key(handles[i]), static_cast<int>(i));
        q.remove(handles[i]);
    }
    ASSERT_EQ(q.size(), values.size() / 2);

    std::sort(values.begin(), values.end(), std::greater<>{});

    for (auto i = 0ul; i < values.size(); i += 2) {
        ASSERT_EQ(q.top(), values[i]);
        q.pop();
    }
    ASSERT_TRUE(q.empty());
}


// in this test, we first fill the priority queue with n elements
// then we change the priority of every second element (to make them larger)
// we pop n/2 elements (i.e. the larger ones) then we reinsert n/2 elements
// and change the keys of all previous elements at last, we pop all n elements
// and see if they are correct (the test shows if it is possible to implement
// workloads that alternate between pushes and pops without compromising handles.
TEST(PriorityQueueTest, CanMixOperations) {
    const std::size_t n = 300;// should be even

    auto q = D::make<std::size_t, 8>(n / 2);
    std::vector<decltype(q)::handle> handles;
    for (std::size_t i = 0; i < n; ++i)
        handles.push_back(q.push(i));
    ASSERT_EQ(q.size(), n);

    for (auto i = 0ul; i < n; i += 2) {
        ASSERT_EQ(q.get_key(handles[i]), i);
        q.change_key(handles[i], 2 * n + i);
    }

    for (auto i = 3 * n - 2; i >= 2 * n; i -= 2) {
        ASSERT_EQ(q.top(), i);
        q.pop();
    }

    ASSERT_EQ(q.size(), n / 2);

    for (auto i = 0ul; i < n; i += 2) {
        handles[i] = q.push(2 * n + i);
    }

    ASSERT_EQ(q.size(), n);

    for (auto i = 1ul; i < n; i += 2) {
        ASSERT_EQ(q.get_key(handles[i]), i);
        q.remove(handles[i]);
    }

    ASSERT_EQ(q.size(), n / 2);

    for (auto i = 3 * n - 2; i >= 2 * n; i -= 2) {
        ASSERT_EQ(q.top(), i);
        q.pop();
    }
    ASSERT_TRUE(q.empty());

    for (std::size_t i = 0; i < n; ++i)
        handles[i] = q.push(i);
    ASSERT_EQ(q.size(), n);

    for (auto i = n; i > 0; --i) {
        ASSERT_EQ(q.top(), i - 1);
        q.pop();
    }
    ASSERT_TRUE(q.empty());
}
