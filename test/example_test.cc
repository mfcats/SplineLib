#include "gmock/gmock.h"

using ::testing::Eq;

TEST(ASimpleTest, IsRunning) {
  ASSERT_THAT(true, Eq(true));
}
