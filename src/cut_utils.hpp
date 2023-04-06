#include <igl/cut_mesh.h>
#include <igl/cut_to_disk.h>
#include <igl/adjacency_list.h>
#include <igl/edges.h>
#include <igl/dijkstra.h>
#include <igl/slim.h>
#include <array>
#include <iostream>
#include <vector>
#include <set>
namespace cut_utils {
    template<typename DerivedV, typename DerivedF>
    void assert_is_disk(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const std::vector<std::vector<int>>& cuts,
        Eigen::PlainObjectBase<DerivedF>&cut_mask) {
        using namespace igl;
        std::set<std::array<int, 2>> cut_edges;
        for (const auto& cut : cuts) {
            const size_t cut_len = cut.size();
            for (size_t i = 0; i < cut_len - 1; i++) {
                std::array<int, 2> e{ cut[i], cut[i + 1] };
                if (e[0] > e[1]) {
                    std::swap(e[0], e[1]);
                }
                cut_edges.insert(e);
            }
        }

        const size_t num_faces = F.rows();
        //Eigen::MatrixXi cut_mask(num_faces, 3);
        //cut_mask.setZero();
        for (size_t i = 0; i < num_faces; i++) {
            std::array<int, 2> e0{ F(i, 0), F(i, 1) };
            std::array<int, 2> e1{ F(i, 1), F(i, 2) };
            std::array<int, 2> e2{ F(i, 2), F(i, 0) };
            if (e0[0] > e0[1]) std::swap(e0[0], e0[1]);
            if (e1[0] > e1[1]) std::swap(e1[0], e1[1]);
            if (e2[0] > e2[1]) std::swap(e2[0], e2[1]);

            if (cut_edges.find(e0) != cut_edges.end()) {
                cut_mask(i, 0) = 1;
            }
            if (cut_edges.find(e1) != cut_edges.end()) {
                cut_mask(i, 1) = 1;
            }
            if (cut_edges.find(e2) != cut_edges.end()) {
                cut_mask(i, 2) = 1;
            }
        }
    }


    void updata_cut_path(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<std::vector<int>>& cuts,
        igl::SLIMData& sData)
    {
        std::vector<int>bound;
        std::set<int>on_bound;
        int nCutNodes = 0;
        int idx_ptr = 0;
        for (auto row : cuts)
        {
            nCutNodes += row.size();
            for (int idx : row)
            {
                bound.push_back(idx);
                on_bound.insert(idx);
            }
        }
        Eigen::VectorXd distortions = sData.E;
        long long int argmax = 0;
        double max_dist = FLT_MIN;
        for (long long int i = 0; i < distortions.size(); i++)
        {

            if (distortions[i] > max_dist)
            {
                
                if (on_bound.find(F(i,0)) != on_bound.end() || 
                    on_bound.find(F(i,1)) != on_bound.end() || 
                    on_bound.find(F(i,2)) != on_bound.end()) continue;
                max_dist = distortions[i];
                argmax = i;
                
            }
        }
        std::vector<std::vector<int>> adj_list;
        igl::adjacency_list(F, adj_list, false);
        int cur_mindis = INT_MAX;
        std::vector<int> shortest_path;
        for (int i = 0; i < bound.size(); i++)
        {
            for (int j = 0; j < F.cols(); j++)
            {
                std::set<int>targets;
                targets.insert(F(argmax, j));

                Eigen::VectorXi min_distance, paths;
                int success = igl::dijkstra(bound[i], targets, adj_list, min_distance, paths);

                if (min_distance(F(argmax, j)) < cur_mindis)
                {
                    shortest_path.clear();
                    igl::dijkstra(F(argmax, j), paths, shortest_path);
                    printf("shortest path length update from %d to %d\n", 
                        cur_mindis, 
                        min_distance(F(argmax, j)));
                    cur_mindis = min_distance(F(argmax, j));
                }

            }
        }
        cuts.push_back(shortest_path);
    }
}