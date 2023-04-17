/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:53:49
 * @LastEditTime: 2023-04-16 17:18:01
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/SubEntity.h
 */
#pragma once
#include <algorithm>
#include <array>

namespace MeshAC{
    class Node;
    class Triangle;
    class Tetrahedron;

    
    class SubEdge{
    public:
        Tetrahedron *tetrahedron=nullptr;
        Triangle *triangle=nullptr;
        int iLocal;
        std::array<Node*, 2> forms;
        std::array<Node*, 2> sortedForms;
        SubEdge(){}
        SubEdge(Node *n0, Node *n1): forms({n0, n1}), sortedForms({n0,n1}){
            std::sort(sortedForms.begin(), sortedForms.end());
        }

    };

    struct SubEdgeHasher
    {
        size_t operator()(const SubEdge &e) const{
            return size_t(e.sortedForms[0]);
        }
    };

    struct SubEdgeEqual
    {
        bool operator()(const SubEdge &e1, const SubEdge &e2) const{
            return (e1.sortedForms[0]==e2.sortedForms[0] && e1.sortedForms[1]==e2.sortedForms[1]);
        }
    };





    class SubTriangle{
    public:
        Tetrahedron *tetrahedron=nullptr;
        int iLocal;
        std::array<Node*, 3> forms;
        std::array<Node*, 3> sortedForms;
        SubTriangle(){};
        SubTriangle(Node *n0, Node *n1, Node *n2): forms({n0, n1, n2}), sortedForms({n0,n1,n2}){
            std::sort(sortedForms.begin(), sortedForms.end());
        }
        bool contained(Node *node){
            return (forms[0]==node || forms[1]==node || forms[2]==node);
        }
    };

    struct SubTriangleHasher
    {
        size_t operator()(const SubTriangle &t) const{
            return size_t(t.sortedForms[0]);
        }
    };

    struct SubTriangleEqual
    {
        bool operator()(const SubTriangle &t1, const SubTriangle &t2) const{
            return (t1.sortedForms[0]==t2.sortedForms[0] 
            && t1.sortedForms[1]==t2.sortedForms[1]
            && t1.sortedForms[2]==t2.sortedForms[2]);
        }
    };
}