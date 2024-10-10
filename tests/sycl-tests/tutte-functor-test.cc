#include <gtest/gtest.h>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/sycl-headers/sycl-parallel-primitives.hh>
#include <thread>
#include <chrono>
#include <future>

/* class DualizeTest : public ::testing::TestWithParam<int> {
protected:
    using T = float;
    using K = uint16_t;
    int N = GetParam();
    Graph G = Graph(neighbours_t(N/2 + 2), true);
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    FullereneBatch<float, uint16_t> batch;
    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
        batch.push_back(G, 0);
        batch.push_back(G, 1);
    }

    void TearDown() override {
        BuckyGen::stop(BQ);
    }
}; */
/* 
class TutteTest : public DualizeTest {};
class SphericalProjectionTest : public DualizeTest {}; */

class FunctorTests : public ::testing::TestWithParam<int> {
protected:
    int N = GetParam();
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    using T = float;
    using K = uint32_t;
    DualizeFunctor<T, K> dualize;
    DualizeFunctor<T, uint16_t> dualize16;

    TutteFunctor<T, K> tutte;
    TutteFunctor<T, uint16_t> tutte16;
    SphericalProjectionFunctor<T, K> spherical_projection;
    SphericalProjectionFunctor<T, uint16_t> spherical_projection16;
    Graph G = Graph(neighbours_t(N/2 + 2), true);
    FullereneBatch<T, K> bigbatch;
    FullereneBatch<T, uint16_t> batch;

    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
        Triangulation C21888(spiral_nomenclature("C21888-[3300,3494,4574,4792,8540,8740,9990,10096,10295,10395,10491,10942]-fullerene"));
        bigbatch.push_back(Graph(C21888.neighbours, true), 0);
        batch.push_back(G, 0);
        batch.push_back(G, 1);
    }

    void TearDown() override {
        BuckyGen::stop(BQ);
    }
};

class MinimumProblem : public ::testing::Test {
protected:
    int N = 20;
    using T = float;
    using K = uint16_t;
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    FullereneBatch<T, K> batch;
    Graph G = Graph(neighbours_t(N/2 + 2), true);
    
    
    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
        batch.push_back(G, 0);
    }

    void TearDown() override {
        BuckyGen::stop(BQ);
    }
};


TEST_P(FunctorTests, AllTestsInOne) {
    auto float_spans_equal = [](auto&& a, auto&& b) {
        return std::equal(a.begin(), a.end(), b.begin(), [](auto a, auto b) {
            auto max_val = std::max(std::abs(a), std::abs(b));
            auto eps = std::numeric_limits<decltype(a)>::epsilon() * 1e3;
            return (std::abs(a - b) / (max_val > eps ? max_val : 1) ) < eps;
        });
    };
    {
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
        });
        ASSERT_EQ(batch1, batch2);
    }
    /* 
    //TODO:: Figure out why dualizing a batch (with small isomers) fails if using a DualizeFunctor with uint32_t as the int type
    {
        FullereneBatch<T, uint32_t> batch1;
        batch1.push_back(G, 0);
        auto batch2 = batch1;
        DualizeFunctor<T, uint32_t> dualize;
        SyclQueue Q(Device::get_devices(DeviceType::CPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
        });
        ASSERT_EQ(batch1, batch2);
    } 
    */
        
    {
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
        });
        ASSERT_TRUE(float_spans_equal(((Span<std::array<float,3>>)(batch1.d_.X_cubic_)).template as_span<float>(), ((Span<std::array<float,3>>)(batch2.d_.X_cubic_)).template as_span<float>() ) );
    }
    {
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SphericalProjectionFunctor<T, uint16_t> spherical_projection;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        spherical_projection(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
            spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
        });
        ASSERT_TRUE(float_spans_equal(((Span<std::array<float,3>>)(batch1.d_.X_cubic_)).template as_span<float>(), ((Span<std::array<float,3>>)(batch2.d_.X_cubic_)).template as_span<float>() ) );
    }

    {
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SphericalProjectionFunctor<T, uint16_t> spherical_projection;
        ForcefieldOptimizeFunctor<PEDERSEN, T, uint16_t> forcefield_optimize;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        spherical_projection(Q, batch1, LaunchPolicy::SYNC);
        forcefield_optimize(Q, batch1, LaunchPolicy::SYNC, 5*N, 5*N);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
            spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
            forcefield_optimize(Q, fullerene, LaunchPolicy::SYNC, 5*N, 5*N);
        });
        ASSERT_TRUE(float_spans_equal(((Span<std::array<float,3>>)(batch1.d_.X_cubic_)).template as_span<float>(), ((Span<std::array<float,3>>)(batch2.d_.X_cubic_)).template as_span<float>() ) );
    }

    {
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SphericalProjectionFunctor<T, uint16_t> spherical_projection;
        ForcefieldOptimizeFunctor<PEDERSEN, T, uint16_t> forcefield_optimize;
        HessianFunctor<PEDERSEN, T, uint16_t> hessian;
        EigenFunctor<EigensolveMode::ENDS, T, uint16_t> eigen;
        EccentricityFunctor<T, uint16_t> eccentricity;
        VolumeFunctor<T, uint16_t> volume;
        SurfaceAreaFunctor<T, uint16_t> surface_area;
        
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        
        SyclVector<T> batch1_eigenvalues(batch.size()*2);
        SyclVector<T> batch1_hessians(batch.size()*N*90);
        SyclVector<uint16_t> batch1_cols(batch.size()*N*90);
        SyclVector<T> batch2_eigenvalues(batch.size()*2);
        SyclVector<T> batch2_hessians(batch.size()*N*90);
        SyclVector<uint16_t> batch2_cols(batch.size()*N*90);
        SyclVector<T> batch1_eigenvectors;
        SyclVector<T> batch1_eccentricities(batch.size());
        SyclVector<T> batch1_volumes(batch.size());
        SyclVector<T> batch1_surface_areas(batch.size());
        SyclVector<T> batch2_eigenvectors;
        SyclVector<T> batch2_eccentricities(batch.size());
        SyclVector<T> batch2_volumes(batch.size());
        SyclVector<T> batch2_surface_areas(batch.size());

        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        spherical_projection(Q, batch1, LaunchPolicy::SYNC);
        forcefield_optimize(Q, batch1, LaunchPolicy::SYNC, 5*N, 5*N);
        hessian(Q, batch1, LaunchPolicy::SYNC, batch1_hessians, batch1_cols);
        eigen(Q, batch1, LaunchPolicy::SYNC, batch1_hessians, batch1_cols, 50, batch1_eigenvalues, batch1_eigenvectors);
        eccentricity(Q, batch1, LaunchPolicy::SYNC, batch1_eccentricities);
        volume(Q, batch1, LaunchPolicy::SYNC, batch1_volumes);
        surface_area(Q, batch1, LaunchPolicy::SYNC, batch1_surface_areas);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
            spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
            forcefield_optimize(Q, fullerene, LaunchPolicy::SYNC, 5*N, 5*N);
            //hessian(Q, fullerene, LaunchPolicy::SYNC, batch2_hessians, batch2_cols);
            //eigen(Q, fullerene, LaunchPolicy::SYNC, batch2_hessians, batch2_cols, 50, batch2_eigenvalues, batch2_eigenvectors);
            eccentricity(Q, fullerene, LaunchPolicy::SYNC, batch2_eccentricities);
            volume(Q, fullerene, LaunchPolicy::SYNC, batch2_volumes);
            surface_area(Q, fullerene, LaunchPolicy::SYNC, batch2_surface_areas);
        });


        auto pol1 = (Polyhedron)batch1[0];
        auto pol2 = (Polyhedron)batch2[0];
        ASSERT_EQ(batch1_eigenvalues.operator Span<float>().subspan(0,2), batch1_eigenvalues.operator Span<float>().subspan(2,2));
        
        ASSERT_FLOAT_EQ(pol1.volume_divergence(), batch1_volumes[0]);
        ASSERT_FLOAT_EQ(pol1.surface_area(), batch1_surface_areas[0]);

        ASSERT_FLOAT_EQ(pol2.volume_divergence(), batch2_volumes[0]);
        ASSERT_FLOAT_EQ(pol2.surface_area(), batch2_surface_areas[0]);
    }
    /* 

    {
        auto batch1 = bigbatch;
        auto batch2 = bigbatch;
        DualizeFunctor<T, K> dualize;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
        });
        Q.wait();
        ASSERT_EQ(batch1, batch2);
    } 

    {
        auto batch1 = bigbatch;
        auto batch2 = bigbatch;
        DualizeFunctor<T, K> dualize;
        TutteFunctor<T, K> tutte;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
        });
        //ASSERT_EQ(batch1, batch2);
        ASSERT_TRUE(float_spans_equal(batch1.d_.X_cubic_, batch2.d_.X_cubic_));
    }
    */


    {
        auto batch1 = bigbatch;
        auto batch2 = bigbatch;
        DualizeFunctor<T, K> dualize;
        TutteFunctor<T, K> tutte;
        SphericalProjectionFunctor<T, K> spherical_projection;
        ForcefieldOptimizeFunctor<PEDERSEN, T, K> forcefield_optimize;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        spherical_projection(Q, batch1, LaunchPolicy::SYNC);
        //forcefield_optimize(Q, batch1, LaunchPolicy::SYNC, 5*N, 5*N);

        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
            spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
            //forcefield_optimize(Q, fullerene, LaunchPolicy::SYNC, 5*N, 5*N);
        });
        ASSERT_TRUE(float_spans_equal(((Span<std::array<float,3>>)(batch1.d_.X_cubic_)).template as_span<float>(), ((Span<std::array<float,3>>)(batch2.d_.X_cubic_)).template as_span<float>() ) );

        //std::cout << batch2[0].d_.A_cubic_.subspan(0, 30) << std::endl;  
        //std::cout << batch1[0].d_.A_cubic_.subspan(0, 30) << std::endl;
        //std::cout << batch2[0].d_.faces_dual_.subspan(0, 30) << std::endl;
        //std::cout << batch1[0].d_.faces_cubic_.subspan(0, 30) << std::endl;
    } 
   
}

/*
TEST_P(FunctorTests, AnotherOne){
        auto batch1 = batch;
        auto batch2 = batch;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        std::for_each(batch2.begin(), batch2.end(), [&](auto fullerene) {
            dualize(Q, fullerene, LaunchPolicy::SYNC);
            tutte(Q, fullerene, LaunchPolicy::SYNC);
        });
        auto float_spans_equal = [](auto& a, auto& b) {
            return std::equal(a.begin(), a.end(), b.begin(), [](auto a, auto b) {
                auto max_val = std::max(std::abs(a), std::abs(b));
                auto eps = std::numeric_limits<decltype(a)>::epsilon() * 1e2;
                return (std::abs(a - b) / (max_val > eps ? max_val : 1) ) < eps;
            });
        };
        ASSERT_TRUE(float_spans_equal(batch1.d_.X_cubic_, batch2.d_.X_cubic_));
}

TEST_P(TutteTest, TutteSmallBatch) {
    //auto batch_copy = batch;
    //DualizeFunctor<T, K> dualize;
    //TutteFunctor<T, K> tutte;
    //SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    //dualize(Q, batch, LaunchPolicy::SYNC);
    //tutte(Q, batch, LaunchPolicy::SYNC);
    std::for_each(batch_copy.begin(), batch_copy.end(), [&](auto fullerene) {
        dualize(Q, fullerene, LaunchPolicy::SYNC);
        tutte(Q, fullerene, LaunchPolicy::SYNC);
    });
    
    ASSERT_TRUE(float_spans_equal(batch.d_.X_cubic_, batch_copy.d_.X_cubic_));
}

TEST_P(SmallIsomerBatchTest, SphericalProjectionSmallBatch) {
    auto batch_copy = batch;
    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    dualize(Q, batch, LaunchPolicy::SYNC);
    tutte(Q, batch, LaunchPolicy::SYNC);
    spherical_projection(Q, batch, LaunchPolicy::SYNC);
    std::for_each(batch_copy.begin(), batch_copy.end(), [&](auto fullerene) {
        dualize(Q, fullerene, LaunchPolicy::SYNC);
        tutte(Q, fullerene, LaunchPolicy::SYNC);
        spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
    });
    ASSERT_EQ(batch_copy, batch);
} */

/* 

TEST_P(FunctorTests, DualizeLargeBatch) {
    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    auto bigbatch_copy = bigbatch;
    dualize(Q, bigbatch, LaunchPolicy::SYNC);

    std::for_each(bigbatch_copy.begin(), bigbatch_copy.end(), [&](auto fullerene) {
        dualize(Q, fullerene, LaunchPolicy::SYNC);
    });

    ASSERT_EQ(bigbatch, bigbatch_copy);
}

TEST_P(FunctorTests, TutteLargeBatch) {
    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    auto bigbatch_copy = bigbatch;
    dualize(Q, bigbatch, LaunchPolicy::SYNC);
    tutte(Q, bigbatch, LaunchPolicy::SYNC);


    std::for_each(bigbatch_copy.begin(), bigbatch_copy.end(), [&](auto fullerene) {
        dualize(Q, fullerene, LaunchPolicy::SYNC);
        tutte(Q, fullerene, LaunchPolicy::SYNC);
    });

    ASSERT_EQ(bigbatch, bigbatch_copy);
}

TEST_P(FunctorTests, SphericalProjectionLargeBatch) {
    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    auto bigbatch_copy = bigbatch;
    dualize(Q, bigbatch, LaunchPolicy::SYNC);
    tutte(Q, bigbatch, LaunchPolicy::SYNC);
    spherical_projection(Q, bigbatch, LaunchPolicy::SYNC);

    std::for_each(bigbatch_copy.begin(), bigbatch_copy.end(), [&](auto fullerene) {
        dualize(Q, fullerene, LaunchPolicy::SYNC);
        tutte(Q, fullerene, LaunchPolicy::SYNC);
        spherical_projection(Q, fullerene, LaunchPolicy::SYNC);
    });

    ASSERT_EQ(bigbatch, bigbatch_copy);
}
*/
INSTANTIATE_TEST_SUITE_P(, FunctorTests, ::testing::Range(20, 21, 20)); 

//INSTANTIATE_TEST_SUITE_P(, DualizeTest, ::testing::Range(20, 21, 20));
//INSTANTIATE_TEST_SUITE_P(, TutteTest, ::testing::Range(20, 21, 20));

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}