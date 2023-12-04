//
// Created by Truong Giang Do on 02/12/2023.
//

#ifndef SSSP_NEW_CUSTOMPRIORITYQUEUE_H
#define SSSP_NEW_CUSTOMPRIORITYQUEUE_H

#include "Graph.h"

template<typename Node>
class custom_priority_queue : public std::priority_queue<Node, std::vector<Node>>
{
public:

bool remove(const Node& value) {
    auto it = std::find(this->c.begin(), this->c.end(), value);

    if (it == this->c.end()) {
        return false;
    }
    if (it == this->c.begin()) {
        // deque the top element
        this->pop();
    }
    else {
        // remove element and re-heap
        this->c.erase(it);
        std::make_heap(this->c.begin(), this->c.end(), this->comp);
    }
    return true;
}
};

#endif //SSSP_NEW_CUSTOMPRIORITYQUEUE_H
