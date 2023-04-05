#pragma once
#include "edge.h"
#include <unordered_map>
#include <vector>

class HashEdgeTable{
public:
    std::unordered_map<int, std::vector<Edge>>  columns;


    HashEdgeTable(){}

    void insert(const Edge &aEdge);
    int key(const Edge &aEdge);
    bool search(Edge &keyEdge, std::vector<Edge> &goalEdge);
    bool exist(const Edge &keyEdge);
    bool remove(int key, const Edge &keyEdge);
    bool remove(const Edge &keyEdge);
};