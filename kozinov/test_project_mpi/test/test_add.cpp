#include "add.h"

#include <gtest.h>


TEST(test_lib, simple_test)
{
  int v  = add(2,2);
  EXPECT_EQ(4, v);
}

