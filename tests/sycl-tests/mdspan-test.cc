#include <gtest/gtest.h>
#include <iostream>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-mdspan.hh>
#include <../../src/sycl/deigen-batch.cc>

TEST(MDSPan, Constructor)
{
  constexpr int m = 6, n = 3, p = 10;
  SyclVector<double> data(m*n*p);
  MDSpan<double, 3> span(data.data(), {m,n,p});
  using array_t = MDSpan<double, 3>::array_t;

  EXPECT_EQ(span.size(), m*n*p);
  EXPECT_EQ(span.shape(), array_t({m,n,p}));
  EXPECT_EQ(span.empty(), false);
}

TEST(MDSPan, StridedAccess)
{
    using array_t = MDSpan<double, 2>::array_t;     
    constexpr int m = 3, n = 2, p = 10;
    SyclVector<double> data(m*p);
    MDSpan<double, 3>  span (data.data(), {m,n,p}, {p*n, p, 1});
    MDSpan<double, 3>  spanT(data.data(), {p,m,n}, {1, p, p*m}); 
    

    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<p;k++){
             span[{i,j,k}] = i*100+j*10+k;
            }
    std::cout << "data = " << data << "\n";

    std::cout << "span.stride = "  << std::span( span.stride().data(),3) << "\n";
    std::cout << "spanT.stride = " << std::span(spanT.stride().data(),3) << "\n";

    
    std::cout << "span = ";
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<p;k++)
                printf("%03g ", span[{i,j,k}]);
    std::cout << "\n";

    std::cout << "spanT = ";
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<p;k++)
                printf("%03g ", spanT[{k,i,j}]);
    std::cout << "\n";
}

TEST(SpanMatrix, CoalesceTest)
{
    constexpr int m = 3, n = 3, p = 4;
    std::array<double,p*m*n> A_data = {1, -1, 1, -1, 1, -1, 1, -1, 0, 0, -1, 0, 1, 1, 1, 1, 1, -1, 1, 0, -1, 0, 1, 0, 1, 1, -1, -1, 0, 0, 0, -1, 1, 1, 1, 0};
    std::array<double,p*m*n> A_coal;
    MDSpan<double,3> A_flat(A_data.data(),{4,3,3}), 
                     A(A_coal.data(), {4,3,3},{1,p*m,p});

    for(int i=0;i<4;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                A[{i,j,k}] = A_flat[{i,j,k}];
    
    for(int i=0;i<4;i++){
        SpanMatrix Ai = A({i,0,0}, 3,3);
        std::cout << "A.shape = " << std::span(Ai.shape().data(),2) << "\n";
        std::cout << "A"<<i << " = \n" << Ai << "\n";
    }
    std::cout << "A.data()" << std::span(A.data(),m*n*p) << "\n";
}

TEST(SpanMatrix, QHQ_test)
{
    constexpr int m = 3, n = 3, p = 4;
    std::array<double,p*m*n> A_data = {0, -1, 2, -1, 2, 1, 2, 1, -2, 0, -1, 0, -1, -2, -1, 0, -1, 2, 0, 0, 1, 0, 2, 0, 1, 0, 0, -2, -2, 1, -2, 0, 1, 1, 1, 2};
    std::array<double,p*m*n> A_coal;
    MDSpan<double,3> A_flat(A_data.data(),{4,3,3}), 
                     A(A_coal.data(), {4,3,3},{1,p*m,p});

    for(int i=0;i<4;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                A[{i,j,k}] = A_flat[{i,j,k}];
    
    double v[n], vc[n], vHA[n];

    for(int i=0;i<4;i++){
        SpanMatrix Ai = A({i,0,0}, 3,3);
        QHQ_workspace w(n, v,vc,vHA);
        std::cout << "Ai = \n" << Ai << "\n";        
        QHQ(Ai, w);
        std::cout << "Ai Tridiagonalized = \n" << Ai << "\n";                
    }
}

#if 0
TEST(SyclVector, ConstructorWithDefaultValue)
{
    SyclVector<int> vec(10, 5);
    EXPECT_EQ(vec.size(), 10);
    EXPECT_EQ(vec.capacity(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
}

TEST(SyclVector, CopyConstructor)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2(vec);
    EXPECT_EQ(vec2.size(), 10);
    EXPECT_EQ(vec2.capacity(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec2[i], 5);
    }
}

TEST(SyclVector, MoveConstructor)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2(std::move(vec));
    EXPECT_EQ(vec2.size(), 10);
    EXPECT_EQ(vec2.capacity(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec2[i], 5);
    }
}

TEST(SyclVector, MoveAssignment)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2;
    vec2 = std::move(vec);
    EXPECT_EQ(vec2.size(), 10);
    EXPECT_EQ(vec2.capacity(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec2[i], 5);
    }
}

TEST(SyclVector, Resize)
{
    SyclVector<int> vec(10, 5);
    vec.resize(20);
    EXPECT_EQ(vec.size(), 20);
    EXPECT_EQ(vec.capacity(), 20);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
}

TEST(SyclVector, ResizeWithDefaultValue)
{
    SyclVector<int> vec(10, 5);
    vec.resize(20, 10);
    EXPECT_EQ(vec.size(), 20);
    EXPECT_EQ(vec.capacity(), 20);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
    for (int i = 10; i < 20; i++)
    {
        EXPECT_EQ(vec[i], 10);
    }
}

TEST(SyclVector, ResizeWithFrontAndBack)
{
    SyclVector<int> vec(10, 5);
    //New size is 20, front is 5, back is 10, seg_size is 5
    vec.resize(20, 5, 10, 5);
    EXPECT_EQ(vec.size(), 20);
    EXPECT_EQ(vec.capacity(), 20);
    for (int i = 0; i < 5; i++)
    {
        //Expect the first 5 elements to be 5
        EXPECT_EQ(vec[i], 5);
    }
    for (int i = 5; i < 20; i++)
    {
        //Expect the rest of the elements to be 0
        EXPECT_EQ(vec[i], 0);
    }
}

TEST(SyclVector, Reserve)
{
    SyclVector<int> vec(10, 5);
    vec.reserve(20);
    EXPECT_EQ(vec.size(), 10);
    EXPECT_EQ(vec.capacity(), 20);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
}

TEST(SyclVector, AssignmentOperator)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2;
    vec2 = vec;
    EXPECT_EQ(vec2.size(), 10);
    EXPECT_EQ(vec2.capacity(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec2[i], 5);
    }
}

TEST(SyclVector, EqualityOperator)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2(10, 5);
    EXPECT_TRUE(vec == vec2);
}

TEST(SyclVector, InequalityOperator)
{
    SyclVector<int> vec(10, 5);
    SyclVector<int> vec2(10, 6);
    EXPECT_FALSE(vec == vec2);
}

TEST(SyclVector, SpanConstructor)
{
    SyclVector<int> vec(10, 5);
    std::span<int> span(vec);
    EXPECT_EQ(span.size(), 10);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(span[i], 5);
    }
}

TEST(SyclVector, PushBack)
{
    SyclVector<int> vec(10, 5);
    vec.push_back(6);
    EXPECT_EQ(vec.size(), 11);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
    EXPECT_EQ(vec[10], 6);
}

TEST(SyclVector, PopBack)
{
    SyclVector<int> vec(10, 5);
    int val = vec.pop_back();
    EXPECT_EQ(val, 5);
    EXPECT_EQ(vec.size(), 9);
    EXPECT_EQ(vec.capacity(), 10);
    for (int i = 0; i < 9; i++)
    {
        EXPECT_EQ(vec[i], 5);
    }
}

TEST(SyclVector, Clear)
{
    SyclVector<int> vec(10, 5);
    vec.clear();
    EXPECT_EQ(vec.size(), 0);
    EXPECT_EQ(vec.capacity(), 10);
}

TEST(SyclVector, BeginEnd)
{
    SyclVector<int> vec(10, 5);
    int i = 0;
    for (auto it = vec.begin(); it != vec.end(); it++)
    {
        EXPECT_EQ(*it, 5);
        i++;
    }
    EXPECT_EQ(i, 10);
}

TEST(SyclVector, At)
{
    SyclVector<int> vec(10, 5);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(vec.at(i), 5);
    }
}
#endif

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
