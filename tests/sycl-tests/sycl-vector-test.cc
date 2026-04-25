#include <gtest/gtest.h>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>


TEST(SyclVector, Constructor)
{
    SyclVector<int> vec(10);
    EXPECT_EQ(vec.size(), 10);
    EXPECT_EQ(vec.capacity(), 10);
}

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

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}