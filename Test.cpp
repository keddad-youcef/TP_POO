#include "gtest/gtest.h"
#include "HypreExample.h"
#include "PETScExample.h"


TEST(UnitTest, MyprojectTestPrintTest){ 
  
  HypreExample HypreExample{};
  PETScExample PETScExample{};

  ResidualNorms ResNormPetObj = PETScExample.run();
  ResidualNorms ResNormHypObj = PETScExample.run();

  // TD3 TESTS
  
  EXPECT_EQ(ResNormPetObj.alien_norm, ResNormPetObj.loc_norm); 
  EXPECT_EQ(ResNormHypObj.alien_norm, ResNormHypObj.loc_norm);


  // TD2 TESTS

  // EXPECT_EQ(PETScExample.run(),0);  
  //EXPECT_NEAR(PETScExample.run() , PETScExample.run(), 0);
}

