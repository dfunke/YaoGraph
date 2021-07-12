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

        int getBalance() const {
            return (left && left->isNode() ? left->asNode()->height : 0) - (right && right->isNode() ? right->asNode()->height : 0);
        }

        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;

        Leaf *leftRep = nullptr;
        Leaf *maxRep = nullptr;
        tIndex height = 0;
    };

    struct Leaf : public Node {

        bool isLeaf() const override { return true; }
        bool isNode() const override { return false; }

        std::unique_ptr<T> obj;

        Leaf *prev = nullptr;
        Leaf *next = nullptr;
    };

public:
    class Iterator {

    protected:
        Leaf *leaf_;

    public:
        Iterator() : leaf_(nullptr) {}
        Iterator(Leaf *leaf) : leaf_(leaf) {}

    public:
        Iterator &operator++() {
            leaf_ = leaf_->next;
            return *this;
        }

        Iterator &operator--() {
            leaf_ = leaf_->prev;
            return *this;
        }

        bool operator==(Iterator b) const { return leaf_ == b.leaf_; }
        bool operator!=(Iterator b) const { return leaf_ != b.leaf_; }

        T &operator*() { return *(leaf_->obj.get()); }
        T *operator->() { return leaf_->obj.get(); }

        Leaf *leaf() { return leaf_; }

        /*operator Leaf *() { return leaf_; }*/
        operator bool() { return leaf_ != nullptr; }

        Iterator prev() {
            return Iterator(leaf_->prev);
        }

        Iterator next() {
            return Iterator(leaf_->next);
        }
    };

public:
    SearchTree() {
        m_root = std::make_unique<InternalNode>();
    }

    Iterator begin() {
        return Iterator(m_first);
    }

    Iterator end() {
        return Iterator(nullptr);
    }

public:
    Iterator insert(Iterator pos, const T &obj) {
        return Iterator(insert(pos.leaf(), obj));
    }

private:
    Leaf *insert(Leaf *pos, const T &obj) {

        // construct new leaf
        std::unique_ptr<Leaf> leaf = std::make_unique<Leaf>();
        leaf->obj = std::make_unique<T>(obj);

        // special case empty list
        if (pos == m_first && pos == nullptr) {

            assert(!m_root->left && !m_root->right);
            assert(m_first == nullptr && m_last == nullptr);

            leaf->prev = nullptr;
            leaf->next = nullptr;
            leaf->parent = m_root.get();

            m_first = leaf.get();
            m_last = leaf.get();

            m_root->left = std::move(leaf);
            updateAndRebalance(m_root.get());

            return m_root->left->asLeaf();
        }

        // take care of the leaf double linked list
        if (pos == end()) {
            leaf->prev = m_last;
            leaf->next = nullptr;

            leaf->prev->next = leaf.get();

            m_last = leaf.get();
        } else if (pos == begin()) {
            leaf->prev = nullptr;
            leaf->next = m_first;

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

            updateAndRebalance(parent);

            return;
        }

        std::unique_ptr<InternalNode> newNode = std::make_unique<InternalNode>();
        newNode->parent = parent;
        newNode->left = std::move(parent->right);
        newNode->left->parent = newNode.get();

        newNode->right = std::move(leaf);
        newNode->right->parent = newNode.get();

        parent->right = std::move(newNode);
        updateAndRebalance(parent->right->asNode());
    }

    void joinFromLeft(InternalNode *parent, std::unique_ptr<Leaf> &&leaf) {

        assert(parent->left);
        if (!parent->right) {
            // parent has empty right child, move left to right, then insert

            leaf->parent = parent;

            parent->right = std::move(parent->left);
            parent->left = std::move(leaf);

            updateAndRebalance(parent);

            return;
        }

        std::unique_ptr<InternalNode> newNode = std::make_unique<InternalNode>();
        newNode->parent = parent;
        newNode->left = std::move(leaf);
        newNode->left->parent = newNode.get();

        newNode->right = std::move(parent->left);
        newNode->right->parent = newNode.get();

        parent->left = std::move(newNode);
        updateAndRebalance(parent->left->asNode());
    }

public:
    Iterator erase(Iterator pos) {
        return Iterator(erase(pos.leaf()));
    }

private:
    Leaf *erase(Leaf *pos) {
        assert(pos != nullptr);
        assert(pos->parent != nullptr);

        // special case singleton list
        if (pos == m_first && pos == m_last) {

            assert(m_root->left && !m_root->right);

            m_first = nullptr;
            m_last = nullptr;

            m_root->left.reset();
            updateAndRebalance(m_root.get());

            return nullptr;
        }

        // take care of the leaf double linked list
        Leaf *retValue = nullptr;// save return value -> node after deleted one

        if (pos == m_last) {
            assert(pos->prev != nullptr);// singleton case already handled above
            pos->prev->next = nullptr;
            m_last = pos->prev;
            retValue = nullptr;
        } else if (pos == m_first) {
            assert(pos->next != nullptr);// singleton case already handled above
            pos->next->prev = nullptr;
            m_first = pos->next;
            retValue = m_first;
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
                parent->right.reset();
            } else {
                assert(child == parent->left.get());
                // we are the left child, move right child into left
                parent->left = std::move(parent->right);
            }

            updateAndRebalance(parent);
        } else {
            assert(parent->left && !parent->right);
            assert(child == parent->left.get());

            // parent becomes empty node
            parent->left.reset();

            if (parent->parent != nullptr) {
                // for non-root delete node recursively
                eraseRec(parent->parent, parent);
            } else {
                // we are root
                assert(parent->isRoot());
                updateAndRebalance(parent);
            }
        }
    }

private:
    // A utility function to right
    // rotate subtree rooted with y
    // See the diagram given above.
    InternalNode *rightRotate(InternalNode *y) {

        std::unique_ptr<Node> yFromParent;
        bool leftChild = false;

        if (y->parent == nullptr) {
            assert(y->isRoot() && y == m_root.get());
            yFromParent = std::move(m_root);
        } else {
            leftChild = (y == y->parent->left.get());
            assert(leftChild || y == y->parent->right.get());
            yFromParent = std::move(leftChild ? y->parent->left : y->parent->right);
        }

        assert(yFromParent && y == yFromParent.get());

        assert(y->left);
        std::unique_ptr<Node> x = std::move(y->left);
        assert(x->isNode() && x->asNode()->left);
        std::unique_ptr<Node> T2 = std::move(x->asNode()->right);

        // perform rotation
        x->parent = y->parent;

        if (T2) {
            y->left = std::move(T2);
            y->left->parent = y;
        } else if (y->right) {
            // if T2 is empty move right child over
            y->left = std::move(y->right);
        } else {
            // y becomes empty, delete it
            yFromParent.reset();
        }

        if (yFromParent) {
            // update info
            updateNodeInfo(y);

            // perform rotation
            x->asNode()->right = std::move(yFromParent);
            x->asNode()->right->parent = x->asNode();
        }

        // update info
        updateNodeInfo(x->asNode());

        // update fromParentPointer
        if (x->parent == nullptr) {
            assert(x->isRoot());
            m_root.reset(static_cast<InternalNode *>(x.release()));

            // Return new root
            return m_root.get();

        } else {
            std::unique_ptr<Node> *pPtr = (leftChild ? &(x->parent->left) : &(x->parent->right));
            (*pPtr) = std::move(x);

            // Return new root
            return (*pPtr)->asNode();
        }
    }

    // A utility function to left
    // rotate subtree rooted with x
    // See the diagram given above.
    InternalNode *leftRotate(InternalNode *x) {

        std::unique_ptr<Node> xFromParent;
        bool leftChild = false;

        if (x->parent == nullptr) {
            assert(x->isRoot() && x == m_root.get());
            xFromParent = std::move(m_root);
        } else {
            leftChild = (x == x->parent->left.get());
            assert(leftChild || x == x->parent->right.get());
            xFromParent = std::move(leftChild ? x->parent->left : x->parent->right);
        }

        assert(xFromParent && x == xFromParent.get());

        assert(x->left && x->right);
        std::unique_ptr<Node> y = std::move(x->right);
        assert(y->isNode() && y->asNode()->left);
        std::unique_ptr<Node> T2 = std::move(y->asNode()->left);

        // perform rotation
        x->right = std::move(T2);
        x->right->parent = x;

        // update info
        updateNodeInfo(x);

        // perform rotation
        y->parent = x->parent;
        y->asNode()->left = std::move(xFromParent);
        y->asNode()->left->parent = y->asNode();

        // update info
        updateNodeInfo(y->asNode());

        // update fromParentPointer
        if (y->parent == nullptr) {
            assert(y->isRoot());
            m_root.reset(static_cast<InternalNode *>(y.release()));

            // Return new root
            return m_root.get();
        } else {
            std::unique_ptr<Node> *pPtr = (leftChild ? &(y->parent->left) : &(y->parent->right));
            (*pPtr) = std::move(y);

            // Return new root
            return (*pPtr)->asNode();
        }
    }

    void updateNodeInfo(InternalNode *node) {
        if (node->left) {
            node->leftRep = (node->left->isNode() ? node->left->asNode()->maxRep : node->left->asLeaf());
            node->maxRep = node->leftRep;
            node->height = 1 + (node->left->isNode() ? node->left->asNode()->height : 0);
        } else {
            assert(!node->right);
            node->leftRep = nullptr;
            node->maxRep = nullptr;
            node->height = 0;
        }

        if (node->right) {
            assert(node->left);
            node->maxRep = (node->right->isNode() ? node->right->asNode()->maxRep : node->right->asLeaf());
            node->height = std::max(node->height, 1 + (node->right->isNode() ? node->right->asNode()->height : 0));
        }
    }

    void updateAndRebalance(InternalNode *node) {

        updateNodeInfo(node);

        int balance = node->getBalance();
        int leftBalance = (node->left && node->left->isNode() ? node->left->asNode()->getBalance() : 0);
        int rightBalance = (node->right && node->right->isNode() ? node->right->asNode()->getBalance() : 0);

        // If this node becomes unbalanced,
        // then there are 4 cases

        // Left Left Case
        if (balance > 1 && leftBalance >= 0) {
            node = rightRotate(node);
        }

        // Left Right Case
        if (balance > 1 && leftBalance < 0) {
            assert(node->left && node->left->isNode());
            leftRotate(node->left->asNode());
            node = rightRotate(node);
        }

        // Right Right Case
        if (balance < -1 && rightBalance <= 0) {
            node = leftRotate(node);
        }

        // Right Left Case
        if (balance < -1 && rightBalance > 0) {
            assert(node->right && node->right->isNode());
            rightRotate(node->right->asNode());
            node = leftRotate(node);
        }


        // recursively traverse up
        if (node->parent != nullptr) {
            updateAndRebalance(node->parent);
        }
    }

public:
    template<class Compare>
    Iterator find(const T &obj, const Compare &cmp) {
        return Iterator(find(m_root.get(), obj, cmp));
    }

private:
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
            return nullptr;
        }
    }

public:
    bool checkInvariants() const {
        return checkInvariants(m_root.get());
    }

private:
    bool checkInvariants(InternalNode *node) const {
        if (std::abs(node->getBalance()) > 1) {
            return false;
        }

        if (!node->left && node->right) {
            return false;
        }

        bool valid = false;
        if (node->left) {
            if (node->left->isNode()) {
                valid = checkInvariants(node->left->asNode());
            } else {
                assert(node->left->isLeaf());
                valid = bool(node->left->asLeaf()->obj);
            }
        }

        if (!valid) { return false; }

        if (node->right) {
            if (node->right->isNode()) {
                valid = checkInvariants(node->right->asNode());
            } else {
                assert(node->right->isLeaf());
                valid = bool(node->right->asLeaf()->obj);
            }
        }

        return valid;
    }

private:
    std::unique_ptr<InternalNode> m_root;

    Leaf *m_first = nullptr;
    Leaf *m_last = nullptr;
};