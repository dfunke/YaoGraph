//
// Created by Daniel Funke on 01.07.21.
//

#pragma once

#include <memory>

#include "Types.hpp"

template<typename T>
class SearchTree {

private:
    // forward declaration of node types
    struct InternalNode;
    struct Leaf;

    // abstract node type
    struct Node {

        virtual ~Node() {}

        bool isRoot() const {
            return parent == nullptr;
        }

        virtual bool isLeaf() const = 0;
        virtual bool isNode() const = 0;

        Leaf *asLeaf() {
            assert(isLeaf());
            return static_cast<Leaf *>(this);
        }

        InternalNode *asNode() {
            assert(isNode());
            return static_cast<InternalNode *>(this);
        }

        // parent is always a internal node
        InternalNode *parent = nullptr;
    };

    struct InternalNode : public Node {

        bool isLeaf() const override { return false; }
        bool isNode() const override { return true; }

        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;

        Leaf *leftRep = nullptr;
        Leaf *maxRep = nullptr;
    };

    struct Leaf : public Node {

        bool isLeaf() const override { return true; }
        bool isNode() const override { return false; }

        std::unique_ptr<T> obj;

        Leaf *prev = nullptr;
        Leaf *next = nullptr;
    };

public:
    SearchTree() {
        m_root = std::make_unique<InternalNode>();
    }

    Leaf *begin() {
        return m_first;
    }

    Leaf *end() {
        return nullptr;
    }

public:
    Leaf *insert(const Leaf *pos, const T &obj) {

        // construct new leaf
        std::unique_ptr<Leaf> leaf = std::make_unique<Leaf>();
        leaf->obj = std::make_unique<T>(obj);

        // special case empty list
        if (pos == begin() && pos == end()) {

            assert(!m_root->left && !m_root->right);
            assert(m_first == nullptr && m_last == nullptr);

            leaf->prev = nullptr;
            leaf->next = end();
            leaf->parent = m_root.get();

            m_first = leaf.get();
            m_last = leaf.get();

            m_root->left = std::move(leaf);
            updateReps(m_root.get());

            return m_root->left->asLeaf();
        }

        // take care of the leaf double linked list
        if (pos == end()) {
            leaf->prev = m_last;
            leaf->next = end();

            leaf->prev->next = leaf.get();

            m_last = leaf.get();
        } else if (pos == begin()) {
            leaf->prev = nullptr;
            leaf->next = begin();

            leaf->next->prev = leaf.get();

            m_first = leaf.get();
        } else {
            leaf->prev = pos->prev;
            assert(leaf->prev != nullptr);
            leaf->next = pos->prev->next;
            assert(leaf->next == pos);

            leaf->prev->next = leaf.get();
            leaf->next->prev = leaf.get();
        }

        // save pointer to leaf, as leaf object looses ownership in join
        Leaf *pLeaf = leaf.get();

        // insert leaf into tree
        if (leaf->prev != nullptr) {
            // we have a left neighbor, join to it from right
            joinFromRight(leaf->prev->parent, std::move(leaf));
        } else {
            assert(leaf->next != nullptr);
            // we have no left neighbor but a right one, join to it from left
            joinFromLeft(leaf->next->parent, std::move(leaf));
        }

        return pLeaf;
    }

private:
    void joinFromRight(InternalNode *parent, std::unique_ptr<Leaf> &&leaf) {

        assert(parent->left);
        if (!parent->right) {
            // parent has empty right child, insert

            leaf->parent = parent;

            parent->right = std::move(leaf);

            updateReps(parent);

            return;
        }

        std::unique_ptr<InternalNode> newNode = std::make_unique<InternalNode>();
        newNode->parent = parent;
        newNode->left = std::move(parent->right);
        newNode->left->parent = newNode.get();

        newNode->right = std::move(leaf);
        newNode->right->parent = newNode.get();

        parent->right = std::move(newNode);
        updateReps(parent->right->asNode());
    }

    void joinFromLeft(InternalNode *parent, std::unique_ptr<Leaf> &&leaf) {

        assert(parent->left);
        if (!parent->right) {
            // parent has empty right child, move left to right, then insert

            leaf->parent = parent;

            parent->right = std::move(parent->left);
            parent->left = std::move(leaf);

            updateReps(parent);

            return;
        }

        std::unique_ptr<InternalNode> newNode = std::make_unique<InternalNode>();
        newNode->parent = parent;
        newNode->left = std::move(leaf);
        newNode->left->parent = newNode.get();

        newNode->right = std::move(parent->left);
        newNode->right->parent = newNode.get();

        parent->left = std::move(newNode);
        updateReps(parent->left->asNode());
    }

public:
    Leaf *erase(Leaf *pos) {
        assert(pos != end());
        assert(pos->parent != nullptr);

        // special case singleton list
        if (pos == m_first && pos == m_last) {

            assert(m_root->left && !m_root->right);

            m_first = nullptr;
            m_last = nullptr;

            m_root->left.release();
            updateReps(m_root.get());

            return end();
        }

        // take care of the leaf double linked list
        Leaf *retValue = nullptr;// save return value -> node after deleted one

        if (pos == m_last) {
            assert(pos->prev != nullptr);// singleton case already handled above
            pos->prev->next = end();
            m_last = pos->prev;
            retValue = end();
        } else if (pos == begin()) {
            assert(pos->next != nullptr);// singleton case already handled above
            pos->next->prev = nullptr;
            m_first = pos->next;
            retValue = begin();
        } else {
            pos->prev->next = pos->next;
            pos->next->prev = pos->prev;
            retValue = pos->next;
        }

        // recursivley delete InternalNodes if they become empty
        eraseRec(pos->parent, pos);

        return retValue;
    }

private:
    void eraseRec(InternalNode *parent, Node *child) {
        assert(parent != nullptr);
        assert(child != nullptr);
        assert(child->parent == parent);

        if (parent->left && parent->right) {
            // parent has both children
            if (child == parent->right.get()) {
                // we are the right child, simply remove
                parent->right.release();
            } else {
                assert(child == parent->left.get());
                // we are the left child, move right child into left
                parent->left = std::move(parent->right);
            }

            updateReps(parent);
        } else {
            assert(parent->left && !parent->right);
            assert(child == parent->left.get());

            // parent becomes empty node
            parent->left.release();

            if (parent->parent != nullptr) {
                // for non-root delete node recursively
                eraseRec(parent->parent, parent);
            } else {
                // we are root
                assert(parent->isRoot());
                updateReps(parent);
            }
        }
    }

private:
    void updateReps(InternalNode *node) {

        if (node->left) {
            node->leftRep = (node->left->isNode() ? node->left->asNode()->maxRep : node->left->asLeaf());
            node->maxRep = node->leftRep;
        } else {
            node->leftRep = nullptr;
            node->maxRep = nullptr;
        }

        if (node->right) {
            assert(node->left);
            node->maxRep = (node->right->isNode() ? node->right->asNode()->maxRep : node->right->asLeaf());
        }

        if (node->parent != nullptr) {
            updateReps(node->parent);
        }
    }

public:
    template<class Compare>
    Leaf *find(const T &obj, const Compare &cmp) {
        return find(m_root.get(), obj, cmp);
    }

    template<class Compare>
    Leaf *find(InternalNode *node, const T &obj, const Compare &cmp) {

        if (cmp(obj, *(node->leftRep->obj))) {
            if (node->left->isNode()) {
                return find(node->left->asNode(), obj, cmp);
            } else {
                assert(node->left->isLeaf());
                return node->left->asLeaf();
            }
        } else if (cmp(obj, *(node->maxRep->obj))) {
            if (node->right->isNode()) {
                return find(node->right->asNode(), obj, cmp);
            } else {
                assert(node->right->isLeaf());
                return node->right->asLeaf();
            }
        } else {
            return end();
        }
    }

private:
    std::unique_ptr<InternalNode> m_root;

    Leaf *m_first = nullptr;
    Leaf *m_last = nullptr;
};