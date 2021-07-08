#include <gtest/gtest.h>

#include <list>

#include "utils/ListIndexTree.hpp"

TEST(ListIndexTreeTest, Adding) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 0 1 2 3 4 5 6 7 8 9 10

    [[maybe_unused]] auto it8 = tree.insert(tree.end(), 8);
    [[maybe_unused]] auto it10 = tree.insert(tree.end(), 10);
    [[maybe_unused]] auto it9 = tree.insert(it10, 9);

    [[maybe_unused]] auto it2 = tree.insert(tree.begin(), 2);
    [[maybe_unused]] auto it1 = tree.insert(tree.begin(), 1);
    [[maybe_unused]] auto it0 = tree.insert(it1, 0);

    for (int i = 3; i < 8; ++i) {
        tree.insert(it8, i);
    }

    // traverse and verify list
    auto it = tree.begin();
    EXPECT_EQ(it, it0);
    EXPECT_EQ(it->prev, nullptr);
    EXPECT_EQ(it->next, it1);
    it = it->next;
    int i = 1;

    for (; it != tree.end(); it = it->next) {
        EXPECT_EQ(it, it->prev->next);

        EXPECT_EQ(*(it->obj), i++);

        if (it->next != nullptr) {
            EXPECT_EQ(it, it->next->prev);
        }
    }
}

TEST(ListIndexTreeTest, Query) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 1 3 5 7 9 11 13 15 17 19
    for (int i = 1; i <= 20; i += 2) {
        tree.insert(tree.end(), i);
    }

    // veriy leaf list
    int v = 1;
    for (auto it = tree.begin(); it != tree.end(); it = it->next) {
        EXPECT_EQ(*(it->obj), v);
        v += 2;
    }

    // traverse and verify list
    for (int i = 0; i < 20; i += 2) {
        auto it = tree.find(i, std::less<int>());

        ASSERT_NE(it, nullptr);
        EXPECT_GT(*(it->obj), i);
        EXPECT_EQ(*(it->obj) - 1, i);

        if (i == 0) {
            EXPECT_EQ(it->prev, nullptr);
        }

        if (i > 0) {
            ASSERT_NE(it->prev, nullptr);
            EXPECT_LT(*(it->prev->obj), i);
        }

        if (i == 18) {
            EXPECT_EQ(it->next, tree.end());
        }
    }

    auto it = tree.find(20, std::less<int>());
    EXPECT_EQ(it, nullptr);
}

TEST(ListIndexTreeTest, Erase) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 0 -- 20
    for (int i = 0; i <= 20; ++i) {
        tree.insert(tree.end(), i);
    }

    // veriy leaf list
    int v = 0;
    for (auto it = tree.begin(); it != tree.end(); it = it->next) {
        EXPECT_EQ(*(it->obj), v++);
    }

    // erase 5- 15; traverse list each time and query for deleted item
    for (int i = 5; i <= 15; ++i) {
        auto it = tree.find(i, std::less<int>());

        ASSERT_NE(it, nullptr);
        ASSERT_NE(it->prev, nullptr);
        EXPECT_GT(*(it->obj), i);
        EXPECT_EQ(*(it->obj), i + 1);
        EXPECT_EQ(*(it->prev->obj), i);

        auto del = tree.erase(it->prev);

        ASSERT_NE(del, nullptr);
        EXPECT_GT(*(del->obj), i);
        EXPECT_EQ(*(del->obj), i + 1);
    }

    v = 0;
    for (auto a = tree.begin(); a != tree.end(); a = a->next) {
        // veriy leaf list
        EXPECT_EQ(*(a->obj), v);

        // check queries
        auto it = tree.find(v, std::less<int>());

        if (v < 20) {
            ASSERT_NE(it, nullptr);
            ASSERT_NE(it->prev, nullptr);
            EXPECT_GT(*(it->obj), v);

            if (v != 4) {
                EXPECT_EQ(*(it->obj), v + 1);
            } else {
                EXPECT_EQ(*(it->obj), 16);
            }

            EXPECT_EQ(*(it->prev->obj), v);
        } else {
            EXPECT_EQ(it, nullptr);
        }


        if (v == 4) {
            v = 16;
        } else {
            ++v;
        }
    }
}