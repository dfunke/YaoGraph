#pragma once

#include <gtest/gtest.h>

constexpr int kInvalid = 0xDEADC0DE;
constexpr int kMovedFrom = 0xBAADBEEF;
constexpr int kValid = 0xCAFEFEED;
constexpr int kEmpty = 0x5EED5AFE;

struct ElementMock {
    static inline int default_constructions = 0;
    static inline int copy_constructions = 0;
    static inline int move_constructions = 0;
    static inline int constructions = 0;
    static inline int copies = 0;
    static inline int moves = 0;
    static inline int destructions = 0;

    static void resetCounters() {
        default_constructions = 0;
        copy_constructions = 0;
        move_constructions = 0;
        constructions = 0;
        copies = 0;
        moves = 0;
        destructions = 0;
    }

    int value = kEmpty;
    int moved = 0;

    ElementMock(int v = kEmpty) : value(v) {
        ++default_constructions;
        ++constructions;
    }

    ElementMock(const ElementMock& rhs) {
        ++copy_constructions;
        ++constructions;
        ++copies;
        value = rhs.value;
    }

    ElementMock(ElementMock&& rhs) noexcept {
        ++move_constructions;
        ++constructions;
        ++moves;
        value = rhs.value;
        moved = rhs.moved + 1;
        rhs.value = kMovedFrom;
    }

    ~ElementMock() {
        ++destructions;
        if (value == kInvalid) ADD_FAILURE() << "Destructor called on invalid element";
        value = kInvalid;
    }

    ElementMock& operator=(const ElementMock& rhs) {
        *this = rhs.value;
        ++copies;
        return *this;
    }

    ElementMock& operator=(ElementMock&& rhs) noexcept {
        *this = rhs.value;
        moved = rhs.moved + 1;
        rhs.value = kMovedFrom;
        ++moves;
        return *this;
    }

    ElementMock& operator=(int v) {
        if (value == kInvalid) ADD_FAILURE() << "Write to invalid element";
        value = v;
        return *this;
    }

    friend bool operator==(const ElementMock& lhs, const ElementMock& rhs) { return lhs.value == rhs.value; }
    friend bool operator!=(const ElementMock& lhs, const ElementMock& rhs) { return lhs.value != rhs.value; }
    friend bool operator<=(const ElementMock& lhs, const ElementMock& rhs) { return lhs.value <= rhs.value; }
    friend bool operator>=(const ElementMock& lhs, const ElementMock& rhs) { return lhs.value >= rhs.value; }
    friend bool operator< (const ElementMock& lhs, const ElementMock& rhs) { return lhs.value <  rhs.value; }
    friend bool operator> (const ElementMock& lhs, const ElementMock& rhs) { return lhs.value >  rhs.value; }

    friend bool operator==(const ElementMock& lhs, int rhs) { return lhs.value == rhs; }
    friend bool operator!=(const ElementMock& lhs, int rhs) { return lhs.value != rhs; }
    friend bool operator<=(const ElementMock& lhs, int rhs) { return lhs.value <= rhs; }
    friend bool operator>=(const ElementMock& lhs, int rhs) { return lhs.value >= rhs; }
    friend bool operator< (const ElementMock& lhs, int rhs) { return lhs.value <  rhs; }
    friend bool operator> (const ElementMock& lhs, int rhs) { return lhs.value >  rhs; }

    friend bool operator==(int lhs, const ElementMock& rhs) { return lhs == rhs.value; }
    friend bool operator!=(int lhs, const ElementMock& rhs) { return lhs != rhs.value; }
    friend bool operator<=(int lhs, const ElementMock& rhs) { return lhs <= rhs.value; }
    friend bool operator>=(int lhs, const ElementMock& rhs) { return lhs >= rhs.value; }
    friend bool operator< (int lhs, const ElementMock& rhs) { return lhs <  rhs.value; }
    friend bool operator> (int lhs, const ElementMock& rhs) { return lhs >  rhs.value; }
};

inline std::ostream& operator<<(std::ostream& os, const ElementMock& v) {
    switch (v.value) {
        case kInvalid: return os << "Invalid";
        case kMovedFrom: return os << "MovedFrom";
        case kValid: return os << "Valid";
        default: return os << v.value;
    }
}
