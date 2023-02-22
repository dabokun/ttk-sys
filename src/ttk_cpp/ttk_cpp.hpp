#include <vector>
#include <string>
class Ttk_rs
{

public:
    Ttk_rs();
    void ttk_helloworld();
    void compute_persistence_pairs(float *data_2d, unsigned int datalen,
                                   unsigned int xdim, unsigned int ydim, int *birth, int *death, int *len,
                                   int *node_ptr, int *arcs_ptr, int *nodes_len, int *arcs_n, int *arcs_len);
    void compute_persistence_pairs_3d(float *data_3d, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim, int *birth, int *death, int *len,
                                      int *node_ptr, int *node_weight_ptr, int *arcs_ptr, int *nodes_len, int *arcs_n, int *arcs_len, int *volume_sizes, int *volume_sizes_len, int thres);
    void data_handling(float *data1, char *data2, float *data3, unsigned int data3_len);
    void get_simplified(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                        float *simplified, int *simplified_len);
    void get_simplified_3d(float *data_dd, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                           float *simplified, int *simplified_len);
    void simplification(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                        unsigned int *cp_point_types, float *cp_coordx, float *cp_coordy, float *cp_value, unsigned int *cp_cellid, unsigned int *cp_pl_vertex_identifier, unsigned int *cp_manifold_size, unsigned int *cp_len,
                        unsigned int *sp_id, float *sp_coordx, float *sp_coordy, unsigned int *sp_point_type, unsigned int *sp_cellid, unsigned int *sp_len,
                        unsigned int *sc_id, unsigned int *sc_source, unsigned int *sc_dest, unsigned int *sc_connectivity_s, unsigned int *sc_connectivity_d, unsigned int *sc_separatrix_id,
                        unsigned int *sc_separatrix_type, unsigned int *sc_f_maxima, unsigned int *sc_f_minima, float *sc_f_diff, unsigned int *sc_len, float *morsesmale);
    void simplification_3d(float *data_3d, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                           unsigned int *cp_point_type, float *cp_coordx, float *cp_coordy, float *cp_coordz, float *cp_value, unsigned int *cp_cellid, unsigned int *cp_pl_vertex_identifier, unsigned int *cp_manifold_size,
                           unsigned int *cp_len, unsigned int *sp_id, float *sp_coordx, float *sp_coordy, float *sp_coordz, unsigned int *sp_point_type, unsigned int *sp_cellid, unsigned int *sp_len,
                           unsigned int *sc_id, unsigned int *sc_source, unsigned int *sc_dest, unsigned int *sc_connectivity_s, unsigned int *sc_connectivity_d, unsigned int *sc_separatrix_id,
                           unsigned int *sc_separatrix_type, unsigned int *sc_f_maxima, unsigned int *sc_f_minima, float *sc_f_diff, unsigned int *sc_len);
    void test_ftr(float *data_3d, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim);
    // void compute_critical_points(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim);
    /*myvec<float, long long int> load_off_file(const std::string &inputPath,
                                              float *pointSet, int pointSetSize,
                                              long long int *trianglesetCo, int trianglesetCoSize,
                                              long long int *trianglesetOff, int trianglesetOffSize);
    void test_ttk_processing1(float *pointSet, int pointSetSize,
                              long long int *trianglesetCo, int trianglesetCoSize,
                              long long int *trianglesetOff, int trianglesetOffSize);
    void test_ttk_processing2(float *pointSet, int pointSetSize);
    int save(const std::vector<float> &pointSet,
             const std::vector<long long int> &triangleSetCo,
             const std::vector<long long int> &triangleSetOff,
             const std::string &outputPath);*/
};
