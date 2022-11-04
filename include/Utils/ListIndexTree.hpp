//
// Created by Daniel Funke on 01.07.21.
//

#pragma once

#include <forward_list>
#include <memory>
#include <vector>

template<typename T>
class SearchTree {

private:
    using index_type = std::size_t;

    template<typename N>
    class MemoryManager {

    public:
        MemoryManager(const index_type &block_size_ = 1e5) : block_size(block_size_) {
            blocks.emplace_front();
            blocks.front().reserve(block_size);
        }

        N *acquire() {
            if (freelist.size()) {
                N *p = freelist[freelist.size() - 1];
                freelist.pop_back();
                return p;
            } else {
                if (blocks.begin()->size() == block_size) {

                    blocks.emplace_front();
                    blocks.front().reserve(block_size);
                }

                blocks.front().emplace_back();
                return &blocks.front()[blocks.front().size() - 1];
            }
        }

        void release(N *p) {
            freelist.push_back(p);
        }

    private:
        std::forward_list<std::vector<N>> blocks;
        index_type block_size;
        std::vector<N *> freelist;
    };

private:
    // forward declaration of node types
    struct InternalNode;
    struct Leaf;

    using height_type = unsigned short;

    // abstract node type
    struct Node {

        virtual ~Node() {}

        bool isRoot() const {
            return parent == nullptr;
        }

        virtual bool isLeaf() const = 0;
        virtual bool isNode() const = 0;

        Leaf *asLeaf() {
            ASSERT(isLeaf());
            return static_cast<Leaf *>(this);
        }

        InternalNode *asNode() {
            ASSERT(isNode());
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

        Node *left;
        Node *right;

        Leaf *leftRep = nullptr;
        Leaf *maxRep = nullptr;
        height_type height = 0;
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

        friend class SearchTree<T>;

    protected:
        Leaf *leaf_;

    public:
        using difference_type = short;
        using value_type = T;
        using pointer = T *;
        using const_pointer = const T *;
        using reference = T &;
        using const_reference = const T &;
        using iterator_category = std::bidirectional_iterator_tag;

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

        reference operator*() { return *(leaf_->obj.get()); }
        pointer operator->() { return leaf_->obj.get(); }

        const_reference operator*() const { return *(leaf_->obj.get()); }
        const_pointer operator->() const { return leaf_->obj.get(); }
    };

public:
    SearchTree() {
        m_root = newInternalNode();

        m_endLeaf = newLeaf();
        m_endLeaf->prev = m_endLeaf;
        m_endLeaf->next = m_endLeaf;

        m_beginLeaf = m_endLeaf;
    }

    Iterator begin() {
        return Iterator(_begin());
    }

    Iterator end() {
        return Iterator(_end());
    }

private:
    Leaf *_begin() {
        return m_beginLeaf;
    }

    Leaf *_end() {
        return m_endLeaf;
    }

public:
    Iterator insert(Iterator pos, const T &obj) {
        return Iterator(insert(pos.leaf_, obj));
    }

private:
    Leaf *newLeaf() {
        return mm_leafs.acquire();
    }

    InternalNode *newInternalNode() {
        return mm_inodes.acquire();
    }

    void deleteNode(Node *p) {
        if (p->isLeaf()) {
            deleteLeaf(p->asLeaf());
        } else {
            deleteInternalNode(p->asNode());
        }
    }

    void deleteLeaf(Leaf *p) {
        mm_leafs.release(p);
    }

    void deleteInternalNode(InternalNode *p) {
        mm_inodes.release(p);
    }

private:
    Leaf *insert(Leaf *pos, const T &obj) {

        // construct new leaf
        Leaf *leaf = newLeaf();
        leaf->obj = std::make_unique<T>(obj);

        // special case empty list
        if (pos == _begin() && pos == _end()) {

            ASSERT(!m_root->left && !m_root->right);
            ASSERT(_begin() == _end());

            leaf->prev = _end();
            leaf->next = _end();
            leaf->parent = m_root;

            m_beginLeaf = leaf;
            m_endLeaf->prev = leaf;

            m_root->left = leaf;
            updateAndRebalance(m_root);

            return m_root->left->asLeaf();
        }

        // take care of the leaf double linked list
        if (pos == _end()) {
            leaf->prev = m_endLeaf->prev;
            leaf->next = _end();

            leaf->prev->next = leaf;

            m_endLeaf->prev = leaf;
        } else if (pos == _begin()) {
            leaf->prev = _end();
            leaf->next = m_beginLeaf;

            leaf->next->prev = leaf;

            m_beginLeaf = leaf;
        } else {
            leaf->prev = pos->prev;
            ASSERT(leaf->prev != _end());
            leaf->next = pos;
            ASSERT(leaf->next == pos);

            leaf->prev->next = leaf;
            leaf->next->prev = leaf;
        }

        // insert leaf into tree
        if (leaf->prev != _end()) {
            // we have a left neighbor, join to it from right
            joinFromRight(leaf->prev->parent, leaf);
        } else {
            ASSERT(leaf->next != _end());
            // we have no left neighbor but a right one, join to it from left
            joinFromLeft(leaf->next->parent, leaf);
        }

        return leaf;
    }

private:
    void joinFromRight(InternalNode *parent, Leaf *leaf) {

        ASSERT(parent->left);
        ASSERT(leaf->prev != _end());

        bool prevIsLeftChild = (parent->left == leaf->prev);
        ASSERT(prevIsLeftChild || (parent->right && parent->right == leaf->prev));

        if (!parent->right) {
            // parent has empty right child, prev must be left child
            // insert new leaf as right child
            ASSERT(prevIsLeftChild);

            leaf->parent = parent;

            parent->right = leaf;

            updateAndRebalance(parent);

            return;
        }

        InternalNode *newNode = newInternalNode();
        newNode->parent = parent;

        newNode->left = prevIsLeftChild ? parent->left : parent->right;
        newNode->left->parent = newNode;

        newNode->right = leaf;
        newNode->right->parent = newNode;

        (prevIsLeftChild ? parent->left : parent->right) = newNode;
        updateAndRebalance((prevIsLeftChild ? parent->left : parent->right)->asNode());
    }

    void joinFromLeft(InternalNode *parent, Leaf *leaf) {

        ASSERT(parent->left);
        // this is only called when object is inserted to the beginning of list
        // new leaf must become leftmost leaf
        ASSERT(leaf->prev == _end());

        if (!parent->right) {
            // parent has empty right child, move left to right, then insert

            leaf->parent = parent;

            parent->right = parent->left;
            parent->left = leaf;

            updateAndRebalance(parent);

            return;
        }

        InternalNode *newNode = newInternalNode();
        newNode->parent = parent;
        newNode->left = leaf;
        newNode->left->parent = newNode;

        newNode->right = parent->left;
        newNode->right->parent = newNode;

        parent->left = newNode;
        updateAndRebalance(parent->left->asNode());
    }

public:
    Iterator erase(Iterator pos) {
        return Iterator(erase(pos.leaf_));
    }

private:
    Leaf *erase(Leaf *pos) {
        ASSERT(pos != nullptr);
        ASSERT(pos != _end());
        ASSERT(pos->parent != nullptr);

        // special case singleton list
        if (pos == _begin() && pos == _end()->prev) {

            ASSERT(m_root->left && !m_root->right);
            ASSERT(m_root->left->isLeaf());

            m_beginLeaf = _end();
            m_endLeaf->prev = _end();

            deleteLeaf(m_root->left->asLeaf());
            m_root->left = nullptr;
            updateAndRebalance(m_root);

            return _end();
        }

        // take care of the leaf double linked list
        Leaf *retValue = nullptr;// save return value -> node after deleted one

        if (pos == _end()->prev) {
            ASSERT(pos->prev != _end());// singleton case already handled above
            pos->prev->next = _end();
            m_endLeaf->prev = pos->prev;
            retValue = _end();
        } else if (pos == _begin()) {
            ASSERT(pos->next != _end());// singleton case already handled above
            pos->next->prev = _end();
            m_beginLeaf = pos->next;
            retValue = _begin();
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
        ASSERT(parent != nullptr);
        ASSERT(child != nullptr);
        ASSERT(child->parent == parent);

        if (parent->left && parent->right) {
            // parent has both children
            if (child == parent->right) {
                // we are the right child, simply remove
                deleteNode(parent->right);
                parent->right = nullptr;
            } else {
                ASSERT(child == parent->left);
                // we are the left child, move right child into left
                parent->left = parent->right;
                parent->right = nullptr;
            }

            ASSERT(parent->left && !parent->right);
            updateAndRebalance(parent);
        } else {
            ASSERT(parent->left && !parent->right);
            ASSERT(child == parent->left);

            // parent becomes empty node
            deleteNode(parent->left);
            parent->left = nullptr;

            if (parent->parent != nullptr) {
                // for non-root delete node recursively
                eraseRec(parent->parent, parent);
            } else {
                // we are root
                ASSERT(parent->isRoot());
                updateAndRebalance(parent);
            }
        }
    }

private:
    // A utility function to right
    // rotate subtree rooted with y
    // See the diagram given above.
    InternalNode *rightRotate(InternalNode *y) {

        //        Node *yFromParent;
        bool leftChild = false;

        if (y->parent == nullptr) {
            ASSERT(y->isRoot() && y == m_root);
            //            yFromParent = m_root;
        } else {
            leftChild = (y == y->parent->left);
            ASSERT(leftChild || y == y->parent->right);
            //            yFromParent = (leftChild ? y->parent->left : y->parent->right);
        }

        //        ASSERT(yFromParent && y == yFromParent);

        ASSERT(y->left);
        Node *x = y->left;
        ASSERT(x->isNode() && x->asNode()->left);
        Node *T2 = x->asNode()->right;

        // perform rotation
        x->parent = y->parent;

        if (T2) {
            y->left = T2;
            y->left->parent = y;
        } else if (y->right) {
            // if T2 is empty move right child over
            y->left = y->right;
            y->right = nullptr;
        } else {
            // y becomes empty, delete it
            deleteInternalNode(y);
            y = nullptr;
        }

        // perform rotation
        x->asNode()->right = y;

        if (y) {
            // update info
            updateNodeInfo(y);
            x->asNode()->right->parent = x->asNode();
        }

        // update info
        updateNodeInfo(x->asNode());

        // update fromParentPointer
        if (x->parent == nullptr) {
            ASSERT(x->isRoot());
            m_root = static_cast<InternalNode *>(x);

            // Return new root
            return m_root;

        } else {
            Node **pPtr = (leftChild ? &(x->parent->left) : &(x->parent->right));
            (*pPtr) = x;

            // Return new root
            return (*pPtr)->asNode();
        }
    }

    // A utility function to left
    // rotate subtree rooted with x
    // See the diagram given above.
    InternalNode *leftRotate(InternalNode *x) {

        //        Node *xFromParent;
        bool leftChild = false;

        if (x->parent == nullptr) {
            ASSERT(x->isRoot() && x == m_root);
            //            xFromParent = m_root;
        } else {
            leftChild = (x == x->parent->left);
            ASSERT(leftChild || x == x->parent->right);
            //            xFromParent = (leftChild ? x->parent->left : x->parent->right);
        }

        //        ASSERT(xFromParent && x == xFromParent);

        ASSERT(x->left && x->right);
        Node *y = x->right;
        ASSERT(y->isNode() && y->asNode()->left);
        Node *T2 = y->asNode()->left;

        // perform rotation
        x->right = T2;
        x->right->parent = x;

        // update info
        updateNodeInfo(x);

        // perform rotation
        y->parent = x->parent;
        y->asNode()->left = x;
        y->asNode()->left->parent = y->asNode();

        // update info
        updateNodeInfo(y->asNode());

        // update fromParentPointer
        if (y->parent == nullptr) {
            ASSERT(y->isRoot());
            m_root = static_cast<InternalNode *>(y);

            // Return new root
            return m_root;
        } else {
            Node **pPtr = (leftChild ? &(y->parent->left) : &(y->parent->right));
            (*pPtr) = y;

            // Return new root
            return (*pPtr)->asNode();
        }
    }

    bool updateNodeInfo(InternalNode *node) {

        auto oldLeftRep = node->leftRep;
        auto oldMaxRep = node->maxRep;
        auto oldHeight = node->height;

        if (node->right) {
            ASSERT(node->left);
            node->leftRep = (node->left->isNode() ? node->left->asNode()->maxRep : node->left->asLeaf());
            node->maxRep = (node->right->isNode() ? node->right->asNode()->maxRep : node->right->asLeaf());
            node->height = 1 + std::max(height_type(node->left->isNode() ? node->left->asNode()->height : 0),
                                        height_type(node->right->isNode() ? node->right->asNode()->height : 0));
        } else if (node->left) {
            ASSERT(!node->right);
            node->leftRep = (node->left->isNode() ? node->left->asNode()->maxRep : node->left->asLeaf());
            node->maxRep = node->leftRep;
            node->height = 1 + (node->left->isNode() ? node->left->asNode()->height : 0);
        } else {
            ASSERT(!node->left && !node->right);
            node->leftRep = nullptr;
            node->maxRep = nullptr;
            node->height = 0;
        }

        return not(oldHeight == node->height && oldLeftRep == node->leftRep && oldMaxRep == node->maxRep);
    }

    void updateAndRebalance(InternalNode *node, bool rebalanceRequired = true) {

        bool nodeInfoChanged = updateNodeInfo(node);

        if (rebalanceRequired) {
            int balance = node->getBalance();

            // If this node becomes unbalanced,
            // then there are 4 cases
            // if this node is unbalanced and we rebalance it, then we do _not_ need to traverese further upward

            if (balance > 1) {
                int leftBalance = (node->left && node->left->isNode() ? node->left->asNode()->getBalance() : 0);

                // Left Left Case
                if (balance > 1 && leftBalance >= 0) {
                    node = rightRotate(node);
                }

                // Left Right Case
                if (balance > 1 && leftBalance < 0) {
                    ASSERT(node->left && node->left->isNode());
                    leftRotate(node->left->asNode());
                    node = rightRotate(node);
                }
            } else if (balance < -1) {
                int rightBalance = (node->right && node->right->isNode() ? node->right->asNode()->getBalance() : 0);

                // Right Right Case
                if (balance < -1 && rightBalance <= 0) {
                    node = leftRotate(node);
                }

                // Right Left Case
                if (balance < -1 && rightBalance > 0) {
                    ASSERT(node->right && node->right->isNode());
                    rightRotate(node->right->asNode());
                    node = leftRotate(node);
                }
            }

            // if we rebalanced the node, we do not need to rebalance further up the tree
            rebalanceRequired = balance <= std::abs(1);
        }

        if (node->parent != nullptr && (nodeInfoChanged || rebalanceRequired)) {
            // recursively traverse up to check whether node above became unbalanced
            updateAndRebalance(node->parent, rebalanceRequired);
        }
    }

public:
    template<typename O, typename Compare>
    Iterator find(const O &obj, const Compare &cmp) {
        return Iterator(find(m_root, obj, cmp));
    }

private:
    template<typename O, typename Compare>
    Leaf *find(InternalNode *node, const O &obj, const Compare &cmp) {

        if (node->leftRep != nullptr) {
            ASSERT(node->maxRep != nullptr);
            ASSERT(node->left);

            if (cmp(obj, *(node->leftRep->obj))) {
                ASSERT(cmp(obj, *(node->maxRep->obj)));
                if (node->left->isNode()) {
                    return find(node->left->asNode(), obj, cmp);
                } else {
                    ASSERT(node->left->isLeaf());
                    ASSERT(cmp(obj, *(node->left->asLeaf()->obj)));
                    return node->left->asLeaf();
                }
            } else if (cmp(obj, *(node->maxRep->obj))) {
                ASSERT(node->right);

                if (node->right->isNode()) {
                    return find(node->right->asNode(), obj, cmp);
                } else {
                    ASSERT(node->right->isLeaf());
                    ASSERT(cmp(obj, *(node->right->asLeaf()->obj)));
                    return node->right->asLeaf();
                }
            } else {
                return _end();
            }
        } else {
            ASSERT(node->isRoot());// only root may have an empty rep when root is empty
            return _end();
        }
    }

public:
    bool checkInvariants() const {
        return checkInvariants(m_root);
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
                ASSERT(node->left->isLeaf());
                valid = bool(node->left->asLeaf()->obj);
            }
        }

        if (!valid) { return false; }

        if (node->right) {
            if (node->right->isNode()) {
                valid = checkInvariants(node->right->asNode());
            } else {
                ASSERT(node->right->isLeaf());
                valid = bool(node->right->asLeaf()->obj);
            }
        }

        return valid;
    }

private:
    InternalNode *m_root;

    Leaf *m_beginLeaf = nullptr;
    Leaf *m_endLeaf;

    MemoryManager<Leaf> mm_leafs;
    MemoryManager<InternalNode> mm_inodes;
};
