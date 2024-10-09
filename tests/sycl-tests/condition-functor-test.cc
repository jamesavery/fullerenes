#include <iostream>
#include <array>
#include <gtest/gtest.h>
#include <fullerenes/sycl-headers/sycl-status-enum.hh>

using namespace condition_detail;

class StatusEnumTest : public ::testing::Test {
protected:
    StatusFlag flag;
};

TEST_F(StatusEnumTest, EmptyFlag)
{
    flag = StatusEnum::EMPTY;
    EXPECT_EQ(flag, StatusEnum::EMPTY);
}

TEST_F(StatusEnumTest, CombineFlags)
{
    flag = StatusFlag::CONVERGED_2D;
    flag |= StatusFlag::FULLERENEGRAPH_PREPARED;
    EXPECT_EQ(flag, (int)StatusEnum::CONVERGED_2D | (int)StatusEnum::FULLERENEGRAPH_PREPARED);
}

TEST_F(StatusEnumTest, TestSingleFlag)
{   
    ConditionFunctor functor(StatusFlag::CONVERGED_2D);
    flag = StatusFlag::CONVERGED_2D;
    EXPECT_EQ(functor(flag), true);
}

TEST_F(StatusEnumTest, TestSingleNotFlag)
{
    ConditionFunctor functor(0);
    flag = 0;
    EXPECT_EQ(functor(flag), true);
}

TEST_F(StatusEnumTest, TestMultipleAndFlags)
{
    ConditionFunctor functor(StatusFlag::CONVERGED_3D | StatusFlag::FULLERENEGRAPH_PREPARED);
    flag = StatusFlag::CONVERGED_3D | StatusFlag::FULLERENEGRAPH_PREPARED;
    EXPECT_EQ(functor(flag), true);
    flag = StatusFlag::DUAL_INITIALIZED;
    EXPECT_EQ(functor(flag), false);
}

TEST_F(StatusEnumTest, TestMultipleOrFlags)
{
    ConditionFunctor functor(0,0, StatusFlag::CONVERGED_2D | StatusFlag::FULLERENEGRAPH_PREPARED);
    flag = StatusFlag::CONVERGED_2D;
    EXPECT_EQ(functor(flag), true);
    flag = StatusFlag::FULLERENEGRAPH_PREPARED | StatusFlag::CONVERGED_2D;
    EXPECT_EQ(functor(flag), true);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
