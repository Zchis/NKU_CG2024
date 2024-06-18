#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"


namespace FluidSimulation
{
    namespace Eulerian3d
    {
        Solver::Solver(MACGrid3d& grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::solve()
        {
            float dt = Eulerian3dPara::dt;

            advect(dt);
            computeforces(dt);
            project(dt);
        }

        void Solver::advect(float dt)
        {
            Glb::GridData3dX newU = mGrid.mU;
            Glb::GridData3dY newV = mGrid.mV;
            Glb::GridData3dZ newW = mGrid.mW;
            Glb::CubicGridData3d newD = mGrid.mD;
            Glb::CubicGridData3d newT = mGrid.mT;

            int numX = Eulerian3dPara::theDim3d[MACGrid3d::X];
            int numY = Eulerian3dPara::theDim3d[MACGrid3d::Y];
            int numZ = Eulerian3dPara::theDim3d[MACGrid3d::Z];

            FOR_EACH_FACE
            {
                if (mGrid.isSolidCell(i, j, k)) {
                    continue;
                }
                glm::vec3 pos_x = mGrid.getBack(i, j, k);
                glm::vec3 new_vel_x = mGrid.semiLagrangian(pos_x, dt);
                newU(i, j, k) = mGrid.getVelocityX(new_vel_x);

                glm::vec3 pos_y = mGrid.getLeft(i, j, k);
                glm::vec3 new_vel_y = mGrid.semiLagrangian(pos_y, dt);
                newV(i, j, k) = mGrid.getVelocityY(new_vel_y);

                glm::vec3 pos_z = mGrid.getBottom(i, j, k);
                glm::vec3 new_vel_z = mGrid.semiLagrangian(pos_z, dt);
                newW(i, j, k) = mGrid.getVelocityZ(new_vel_z);
                
            }

            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j, k)) {
                    continue;
                }
                glm::vec3 pos_p = mGrid.getCenter(i, j, k);
                glm::vec3 new_vel_p = mGrid.semiLagrangian(pos_p, dt);

                newD(i, j, k) = mGrid.getDensity(new_vel_p);
                newT(i, j, k) = mGrid.getTemperature(new_vel_p);
            }

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mW = newW;
            mGrid.mD = newD;
            mGrid.mT = newT;
        }

        void Solver::computeforces(float dt)
        {
            int numX = mGrid.dim[0];
            int numY = mGrid.dim[1];
            int numZ = mGrid.dim[2];

            Glb::GridData3dZ newW = mGrid.mW;
            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j, k) || mGrid.isSolidCell(i, j, k + 1) || mGrid.isSolidCell(i, j, k - 1)) {
                    continue;
                }

                glm::vec3 pos1 = mGrid.getCenter(i, j, k);
                glm::vec3 pos0 = mGrid.getCenter(i, j, k - 1);

                float bforce1 = mGrid.getBoussinesqForce(pos1);
                float bforce0 = mGrid.getBoussinesqForce(pos0);

                float w = (bforce1 + bforce0) * 0.5 * dt;
                newW(i, j, k) += w;
            }

            mGrid.mW = newW;
        }

        void Solver::project(float dt)
        {
            int numX = mGrid.dim[0];
            int numY = mGrid.dim[1];
            int numZ = mGrid.dim[2];

            Glb::GridData3d newP = mGrid.mD;
            newP.initialize(0.0);
            Glb::GridData3dY newV = mGrid.mV;
            Glb::GridData3dX newU = mGrid.mU;
            Glb::GridData3dZ newW = mGrid.mW;

            float airDensity = Eulerian3dPara::airDensity;
            float cellSize = mGrid.cellSize;


            /*
            for (int iteratenum = 40; iteratenum > 0; iteratenum--) {
                FOR_EACH_FACE
                {
                    if (mGrid.isSolidCell(i, j, k)) {
                        continue;
                    }
                    double s = mGrid.getPressureCoeffBetweenCells(i, j, k, i, j, k);
                    if (s == 0) {
                        continue;
                    }
                    double cp = airDensity * cellSize / dt;
                    double div = (newU(i + 1, j, k) - newU(i, j, k) + newV(i, j + 1, k) - newV(i, j, k) + newW(i, j, k + 1) - newW(i, j, k));
                    double p = -div / s;

                    newU(i, j, k) -= p * (!mGrid.mSolid(i - 1, j, k));
                    newU(i + 1, j, k) += p * (!mGrid.mSolid(i + 1, j, k));
                    newV(i, j, k) -= p * (!mGrid.mSolid(i, j - 1, k));
                    newV(i, j + 1, k) += p * (!mGrid.mSolid(i, j + 1, k));
                    newW(i, j, k) -= p * (!mGrid.mSolid(i, j, k - 1));
                    newW(i, j, k + 1) += p * (!mGrid.mSolid(i, j, k + 1));
                }
            }

            for (int i = 0; i <= numX; i++) {
                for (int k = 0; k <= numZ; k++) {
                    newV(i, numY - 1, k) = 0;
                }
            }
            for (int j = 0; j <= numY; j++) {
                for (int k = 0; k <= numZ; k++) {
                    newU(numX - 1, j, k) = 0;
                } 
            }
            for (int i = 0; i <= numX; i++) {
                for (int j = 0; j <= numY; j++) {
                    newW(i, j, numZ - 1) = 0;
                }
            }

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mW = newW;
            */

                // 迭代计算压力
            for (int iteration = 100; iteration > 0; iteration--) {
                FOR_EACH_CELL{
                    if (mGrid.isSolidCell(i, j, k)) {
                        continue;
                    }

                    double px1 = mGrid.isSolidCell(i + 1, j, k) ? 0.0 : newP(i + 1, j, k);
                    double px0 = mGrid.isSolidCell(i - 1, j, k) ? 0.0 : newP(i - 1, j, k);

                    double py1 = mGrid.isSolidCell(i, j + 1, k) ? 0.0 : newP(i, j + 1, k);
                    double py0 = mGrid.isSolidCell(i, j - 1, k) ? 0.0 : newP(i, j - 1, k);

                    double pz1 = mGrid.isSolidCell(i, j, k + 1) ? 0.0 : newP(i, j, k + 1);
                    double pz0 = mGrid.isSolidCell(i, j, k - 1) ? 0.0 : newP(i, j, k - 1);

                    double div = mGrid.getDivergence(i, j, k);

                    double b = -1 * div * airDensity * cellSize * cellSize / dt;
                    double sum = px1 + px0 + py1 + py0 + pz1 + pz0;
                    double s = mGrid.getPressureCoeffBetweenCells(i, j, k, i, j, k);
                    newP(i, j, k) = (b + sum) / s;
                }
            }

            // 调整速度场
            FOR_EACH_CELL{
                newU(i + 1, j, k) -= dt * (newP(i + 1, j, k) - newP(i, j, k)) / (cellSize * airDensity);
                newV(i, j + 1, k) -= dt * (newP(i, j + 1, k) - newP(i, j, k)) / (cellSize * airDensity);
                newW(i, j, k + 1) -= dt * (newP(i, j, k + 1) - newP(i, j, k)) / (cellSize * airDensity);
            }

            // 边界处理
            FOR_EACH_CELL
            {
                if (mGrid.isSolidFace(i, j, k, MACGrid3d::Direction::X)) {
                    newU(i, j, k) = 0;
                }
                if (mGrid.isSolidFace(i, j, k, MACGrid3d::Direction::Y)) {
                    newV(i, j, k) = 0;
                }
                if (mGrid.isSolidFace(i, j, k, MACGrid3d::Direction::Z)) {
                    newW(i, j, k) = 0;
                }
            }

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mW = newW;

        }
    }
}
