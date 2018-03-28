#include "gmock/gmock.h"

using ::testing::Eq;

TEST(ASimpleTest, IsRunning) { // NOLINT
  ASSERT_THAT(true, Eq(true));
}
