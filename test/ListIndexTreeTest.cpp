#include <gtest/gtest.h>

#include <list>

#include "utils/ASSERT.hpp"
#include "utils/ListIndexTree.hpp"

TEST(ListIndexTreeTest, Adding) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 0 1 2 3 4 5 6 7 8 9 10

    [[maybe_unused]] auto it8 = tree.insert(tree.end(), 8);
    EXPECT_TRUE(tree.checkInvariants());
    [[maybe_unused]] auto it10 = tree.insert(tree.end(), 10);
    EXPECT_TRUE(tree.checkInvariants());
    [[maybe_unused]] auto it9 = tree.insert(it10, 9);
    EXPECT_TRUE(tree.checkInvariants());

    [[maybe_unused]] auto it2 = tree.insert(tree.begin(), 2);
    EXPECT_TRUE(tree.checkInvariants());
    [[maybe_unused]] auto it1 = tree.insert(tree.begin(), 1);
    EXPECT_TRUE(tree.checkInvariants());
    [[maybe_unused]] auto it0 = tree.insert(it1, 0);
    EXPECT_TRUE(tree.checkInvariants());

    for (int i = 3; i < 8; ++i) {
        tree.insert(it8, i);
        EXPECT_TRUE(tree.checkInvariants());
    }

    // traverse and verify list
    auto it = tree.begin();
    EXPECT_EQ(it, it0);
    EXPECT_EQ(std::prev(it), tree.end());
    EXPECT_EQ(std::next(it), it1);
    it = std::next(it);
    int i = 1;

    for (; it != tree.end(); ++it) {
        EXPECT_EQ(it, std::next(std::prev(it)));

        EXPECT_EQ(*it, i++);

        EXPECT_EQ(it, std::prev(std::next(it)));
    }
}

TEST(ListIndexTreeTest, Query) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 1 3 5 7 9 11 13 15 17 19
    for (int i = 1; i <= 20; i += 2) {
        tree.insert(tree.end(), i);
        EXPECT_TRUE(tree.checkInvariants());
    }

    // veriy leaf list
    int v = 1;
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        EXPECT_EQ(*it, v);
        v += 2;
    }

    // traverse and verify list
    for (int i = 0; i < 20; i += 2) {
        auto it = tree.find(i, std::less<int>());

        ASSERT_NE(it, tree.end());
        EXPECT_GT(*it, i);
        EXPECT_EQ(*it - 1, i);

        if (i == 0) {
            EXPECT_EQ(std::prev(it), tree.end());
        }

        if (i > 0) {
            ASSERT_NE(std::prev(it), tree.end());
            EXPECT_LT(*std::prev(it), i);
        }

        if (i == 18) {
            EXPECT_EQ(std::next(it), tree.end());
        }
    }

    auto it = tree.find(20, std::less<int>());
    EXPECT_EQ(it, tree.end());
}

TEST(ListIndexTreeTest, Erase) {

    SearchTree<int> tree;

    // build tree that should have leaf list like: 0 -- 20
    for (int i = 0; i <= 20; ++i) {
        tree.insert(tree.end(), i);
        EXPECT_TRUE(tree.checkInvariants());
    }

    // veriy leaf list
    int v = 0;
    for (auto it = tree.begin(); it != tree.end(); it = ++it) {
        EXPECT_EQ(*it, v++);
    }

    // erase 5- 15; traverse list each time and query for deleted item
    for (int i = 5; i <= 15; ++i) {
        auto it = tree.find(i, std::less<int>());

        ASSERT_NE(it, tree.end());
        ASSERT_NE(std::prev(it), tree.end());
        EXPECT_GT(*it, i);
        EXPECT_EQ(*it, i + 1);
        EXPECT_EQ(*std::prev(it), i);

        auto del = tree.erase(std::prev(it));
        EXPECT_TRUE(tree.checkInvariants());

        ASSERT_NE(del, tree.end());
        EXPECT_GT(*del, i);
        EXPECT_EQ(*del, i + 1);
    }

    v = 0;
    for (auto a = tree.begin(); a != tree.end(); ++a) {
        // veriy leaf list
        EXPECT_EQ(*a, v);

        // check queries
        auto it = tree.find(v, std::less<int>());

        if (v < 20) {
            ASSERT_NE(it, tree.end());
            ASSERT_NE(std::prev(it), tree.end());
            EXPECT_GT(*it, v);

            if (v != 4) {
                EXPECT_EQ(*it, v + 1);
            } else {
                EXPECT_EQ(*it, 16);
            }

            EXPECT_EQ(*std::prev(it), v);
        } else {
            EXPECT_EQ(it, tree.end());
        }


        if (v == 4) {
            v = 16;
        } else {
            ++v;
        }
    }
}