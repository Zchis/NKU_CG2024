#include "Eulerian/include/Solver.h"
#include "Configure.h"
#include <omp.h>

/*
namespace FluidSimulation
{
    namespace Eulerian2d
    {
        Solver::Solver(MACGrid2d& grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::solve()
        {
            // TODO
            // Solves the fluid simulation by performing some steps, which may include:
            // 1. advection
            // 2. compute external forces
            // 3. projection
            // ...

        }
    }
}
*/

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        Solver::Solver(MACGrid2d& grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::solve()
        {
            float dt = Eulerian2dPara::dt;

            // ��һ��: ����
            advect(dt);

            // �ڶ���: �����ⲿ��
            computeforces(dt);

            // ������: ͶӰ
            project(dt);

        }

        void Solver::advect(float dt)
        {
            // �����������P
            // ʹ�ð��������շ���
            Glb::GridData2dX newU = mGrid.mU;
            Glb::GridData2dY newV = mGrid.mV;
            Glb::CubicGridData2d newD = mGrid.mD;
            Glb::CubicGridData2d newT = mGrid.mT;

            int numX = Eulerian2dPara::theDim2d[MACGrid2d::X];
            int numY = Eulerian2dPara::theDim2d[MACGrid2d::Y];

            // �����ٶ�
            FOR_EACH_LINE{
                // �ж��ǹ�����߽߱�
                if (mGrid.isSolidCell(i, j)) {
                    continue;
                }
                
                // V
                if (i < numX) {
                    glm::vec2 pos_y = mGrid.getBottom(i, j);
                    glm::vec2 new_vel_y = mGrid.semiLagrangian(pos_y, dt);
                    // glm::vec2 new_vel_y = mGrid.RK2(pos_y, dt);
                    newV(i, j) = mGrid.getVelocityY(new_vel_y);
                }
                
                // U
                
                if (j < numY) {
                    glm::vec2 pos_x = mGrid.getLeft(i, j);
                    glm::vec2 new_vel_x = mGrid.semiLagrangian(pos_x, dt);
                    // glm::vec2 new_vel_x = mGrid.RK2(pos_x, dt);
                    newU(i, j) = mGrid.getVelocityX(new_vel_x);
                }
                
            }

            // ��������
            FOR_EACH_CELL
            {
                // �ж��ǹ�����߽߱�
                if (mGrid.isSolidCell(i, j)) {
                    continue;
                }
                glm::vec2 pos_p = mGrid.getCenter(i, j);
                glm::vec2 new_vel_p = mGrid.semiLagrangian(pos_p, dt);
                // glm::vec2 new_vel_p = mGrid.RK2(pos_p, dt);
                
                newD(i, j) = mGrid.getDensity(new_vel_p);
                newT(i, j) = mGrid.getTemperature(new_vel_p);
            }

            // �߽�����

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mD = newD;
            mGrid.mT = newT;
        }

        void Solver::computeforces(float dt)
        {
            int numX = mGrid.dim[0];
            int numY = mGrid.dim[1];
            // ����
            Glb::GridData2dY newV = mGrid.mV;
            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j) || mGrid.isSolidCell(i, j - 1) || mGrid.isSolidCell(i, j + 1)) {
                    continue;
                }
                
                glm::vec2 pos1 = mGrid.getCenter(i, j);
                glm::vec2 pos0 = mGrid.getCenter(i, j - 1);
                // ���ϵ�������
                float bforce1 = mGrid.getBoussinesqForce(pos1);
                float bforce0 = mGrid.getBoussinesqForce(pos0);

                float v = (bforce1 + bforce0) * 0.5 * dt;
                // ���� v ����
                newV(i, j) += v;
            }

            mGrid.mV = newV;

        }

        void Solver::project(float dt)
        {
            int numX = mGrid.dim[0];
            int numY = mGrid.dim[1];
            Glb::CubicGridData2d newP = mGrid.mP;
            newP.initialize(0.0);
            Glb::GridData2dY newV = mGrid.mV;
            Glb::GridData2dX newU = mGrid.mU;

            float aird = Eulerian2dPara::airDensity;

            float cellSize = mGrid.cellSize;

            // first
            
            /*
            double cp = aird * cellSize / dt;
            for (int iteratenum = 40; iteratenum > 0; iteratenum--) {
                FOR_EACH_CELL
                {
                    
                    if (mGrid.isSolidCell(i, j)) {
                        continue;
                    }
                    double s = mGrid.getPressureCoeffBetweenCells(i, j, i, j);
                    if (s == 0) {
                        continue;
                    }
                    double div = (newU(i + 1, j) - newU(i, j) + newV(i, j + 1) - newV(i, j));
                    double p = - div / s;

                    newP(i, j) += cp * p;

                    newU(i, j) -= p * (!mGrid.mSolid(i - 1, j));
                    newU(i + 1, j) += p * (!mGrid.mSolid(i + 1, j));
                    newV(i, j) -= p * (!mGrid.mSolid(i, j - 1));
                    newV(i, j + 1) += p * (!mGrid.mSolid(i, j + 1));
                
                };
            }
            
            FOR_EACH_LINE
            {
                // ��U
                if (mGrid.isSolidFace(i, j, MACGrid2d::Direction::X)) {
                    newU(i, j) = 0;
                }
                // ��V
                if (mGrid.isSolidFace(i, j, MACGrid2d::Direction::Y)) {
                    newV(i, j) = 0;
                }
            }
            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mP = newP;         
            */

            // second
            for (int iteration = 100; iteration > 0; iteration--) {
                FOR_EACH_CELL{
                    if (mGrid.isSolidCell(i, j)) {
                        continue;
                    }
                    /*
                    if (mGrid.isSolidCell(i - 1, j)) {
                        newP(i - 1, j) = newP(i, j) - cellSize * aird * newU(i + 1, j) / dt;
                    }
                    if (mGrid.isSolidCell(i, j - 1)) {
                        newP(i, j - 1) = newP(i, j) - cellSize * aird * newV(i, j + 1) / dt;
                    }
                    */ 
                    double px1 = mGrid.isSolidCell(i + 1, j) ? 0.0 : newP(i + 1, j);
                    double px0 = mGrid.isSolidCell(i - 1, j) ? 0.0 : newP(i - 1, j);

                    double py1 = mGrid.isSolidCell(i, j + 1) ? 0.0 : newP(i, j + 1);
                    double py0 = mGrid.isSolidCell(i, j - 1) ? 0.0 : newP(i, j - 1);
                    

                    double div = mGrid.getDivergence(i, j);

                    // b
                    // double b = -1 * (newU(i + 1, j) - newU(i, j) + newV(i, j + 1) - newV(i, j)) * (aird) * cellSize / (dt);
                    double b = -1 * (div) * (aird) * cellSize * cellSize / (dt);
                    // sum
                    double sum = (px1 + px0 + py1 + py0);
                    double s = mGrid.getPressureCoeffBetweenCells(i, j, i, j);
                    newP(i, j) = (b + sum) / s;
                };
            }
           
            FOR_EACH_CELL{
                newU(i + 1,j) -= dt * (newP(i + 1,j) - newP(i,j)) / (cellSize * aird);
                newV(i, j + 1) -= dt * (newP(i, j + 1) - newP(i, j)) / (cellSize * aird);
            }

            // �߽紦��
            FOR_EACH_LINE
            {
                // ��U
                if (mGrid.isSolidFace(i, j, MACGrid2d::Direction::X)) {
                    newU(i, j) = 0;
                }
                // ��V
                if (mGrid.isSolidFace(i, j, MACGrid2d::Direction::Y)) {
                    newV(i, j) = 0;
                }
            }
            
            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mP = newP;


            // third
            // �����������ճ��ϵ��
            // ���ȼ����м��ٶ�
            // Ԥ�㲽 Step1��
            /*
            FOR_EACH_CELL{
                if (mGrid.isSolidCell(i, j)) {
                    continue;
                }
                if (mGrid.isSolidCell(i - 1, j) == 1) {
                    newU(i - 1, j) = mGrid.mU(i, j);
                }
                if (mGrid.isSolidCell(i, j - 1) == 1) {
                    newV(i, j - 1) = mGrid.mV(i, j);
                }
                //newU
                double ududx = mGrid.mU(i,j) * (mGrid.mU(i + 1,j) - mGrid.mU(i - 1,j)) / (2 * cellSize);
                double vdudy = 0.25 * (mGrid.mV(i - 1, j) + mGrid.mV(i, j) + mGrid.mV(i - 1, j + 1) + mGrid.mV(i, j + 1)) * (mGrid.mU(i, j + 1) - mGrid.mU(i, j - 1)) / (2 * cellSize);
                newU(i, j) = mGrid.mU(i, j) + dt * (-(ududx + vdudy));
                //newV
                double udvdx = 0.25 * (mGrid.mU(i, j) + mGrid.mU(i, j - 1) + mGrid.mU(i + 1, j) + mGrid.mU(i + 1, j - 1)) * (mGrid.mV(i + 1, j) - mGrid.mV(i - 1, j)) / (2 * cellSize);
                double vdvdy = mGrid.mV(i, j) * (mGrid.mV(i, j + 1) - mGrid.mV(i, j - 1)) / (2 * cellSize);
                newV(i, j) = mGrid.mV(i, j) + dt * (-(udvdx + vdvdy));
            };

            // ѹ��Poisson ����
            // ������ܶȻ�仯
            for (int iteration = 0; iteration < 100; iteration++) {
                FOR_EACH_CELL{
                    if (mGrid.isSolidCell(i, j)) {
                        continue;
                    }
                    if (mGrid.isSolidCell(i - 1, j) == 1) {
                        newP(i - 1, j) = newP(i, j);
                    }
                    if (mGrid.isSolidCell(i, j - 1) == 1) {
                        newP(i, j - 1) = newP(i, j);
                    }
                    // Gauss-Seidel ��������ѹ��
                    double dudx = (newU(i + 1,j) - newU(i,j)) / cellSize;
                    double dvdy = (newV(i, j + 1) - newV(i, j)) / cellSize;
                    double R = (dudx + dvdy) * mGrid.mD(i, j) / dt;

                    double p1 = newP(i - 1, j) / (cellSize * cellSize);
                    double p2 = newP(i + 1, j) / (cellSize * cellSize);
                    double p3 = newP(i, j - 1) / (cellSize * cellSize);
                    double p4 = newP(i, j + 1) / (cellSize * cellSize);

                    newP(i, j) = (p1 + p2 + p3 + p4 - R) / (2 / (cellSize * cellSize) + 2 / (cellSize * cellSize));

                }
            }

            // ������õ�ѹ�����м��ٶ������һʱ�̵��ٶȳ�
            FOR_EACH_CELL{
                if (mGrid.isSolidCell(i, j)) {
                    continue;
                }
                if (mGrid.isSolidCell(i - 1, j) == 1) {
                    newU(i - 1, j) = mGrid.mU(i, j);
                }
                if (mGrid.isSolidCell(i, j - 1) == 1) {
                    newV(i, j - 1) = mGrid.mV(i, j);
                }
                newU(i, j) = newU(i, j) - (dt / aird) * (newP(i, j) - newP(i - 1, j)) / cellSize;
                newV(i, j) = newV(i, j) - (dt / aird) * (newP(i, j) - newP(i, j - 1)) / cellSize;
            }

            // �߽���������
            // �߽紦��
            setbnd(1, numX, numY, newU);
            setbnd(2, numX, numY, newV);

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mP = newP;
            
            */
            
            // forth
            /*
            for (int iteratenum = 40; iteratenum > 0; iteratenum--)
            {
                // ����ż����ż���к�������������
            #pragma omp parallel for collapse(2)
                for (int j = 0; j < numY; j++)
                {
                    for (int i = 0; i < numX; i++)
                    {
                        if ((i + j) % 2 == 0)  // ż����ż���� �� ������������
                        {
                            if (mGrid.isSolidCell(i, j))
                            {
                                continue;
                            }

                            double s = mGrid.getPressureCoeffBetweenCells(i, j, i, j);
                            if (s == 0)
                            {
                                continue;
                            }
                            double cp = airDensity * cellSize / dt;
                            double div = (newU(i + 1, j) - newU(i, j) + newV(i, j + 1) - newV(i, j));
                            double p = -div / s;

                            newP(i, j) += cp * p;

                            newU(i, j) -= p * (!mGrid.mSolid(i - 1, j));
                            newU(i + 1, j) += p * (!mGrid.mSolid(i + 1, j));
                            newV(i, j) -= p * (!mGrid.mSolid(i, j - 1));
                            newV(i, j + 1) += p * (!mGrid.mSolid(i, j + 1));
                        }
                    }
                }

                // ����������ż���к�ż����������
            #pragma omp parallel for collapse(2)
                for (int j = 0; j < numY; j++)
                {
                    for (int i = 0; i < numX; i++)
                    {
                        if ((i + j) % 2 != 0)  // ������ż���� �� ż����������
                        {
                            if (mGrid.isSolidCell(i, j))
                            {
                                continue;
                            }

                            double s = mGrid.getPressureCoeffBetweenCells(i, j, i, j);
                            if (s == 0)
                            {
                                continue;
                            }
                            double cp = airDensity * cellSize / dt;
                            double div = (newU(i + 1, j) - newU(i, j) + newV(i, j + 1) - newV(i, j));
                            double p = -div / s;

                            newP(i, j) += cp * p;

                            newU(i, j) -= p * (!mGrid.mSolid(i - 1, j));
                            newU(i + 1, j) += p * (!mGrid.mSolid(i + 1, j));
                            newV(i, j) -= p * (!mGrid.mSolid(i, j - 1));
                            newV(i, j + 1) += p * (!mGrid.mSolid(i, j + 1));
                        }
                    }
                }
            }

            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mP = newP;
            */
        }

        /*
        void Solver::setbnd(int type, int numX, int numY, Glb::GridData2d& data) {
            // �߽���������
            // �߽紦��
            
            // type == 0: ���P
            // type == 1: ��U
            // type == 2: ��V
            if (type == 0) {
                for (int i = 0; i < numX; i++) {
                    data(i, -1) = data(i, 0);
                    data(i, numY) = data(i, numY - 1);
                }
                for (int j = 0; j < numY; j++) {
                    data(-1, j) = data(0, j);
                    data(numX, j) = data(numX - 1, j);
                }
            }
            if (type == 1) {
                for (int i = 0; i < numX; i++) {
                    data(i, -1) = data(i, 0);
                    data(i, numY) = data(i, numY - 1);
                }
                for (int j = 0; j < numY; j++) {
                    data(-1, j) = -data(0, j);
                    data(numX, j) = -data(numX - 1, j);
                }
            }
            if (type == 2) {
                for (int i = 0; i < numX; i++) {
                    data(i, -1) = -data(i, 0);
                    data(i, numY) = -data(i, numY - 1);           
                }
                for (int j = 0; j < numY; j++) {
                    data(-1, j) = data(0, j);
                    data(numX, j) = data(numX - 1, j);
                }
            }
            data(-1, -1) = 0.5 * (data(0, -1) + data(-1, 0));
            data(-1, numY) = 0.5 * (data(-1, numY - 1) + data(0, numY - 1));
            data(numX, -1) = 0.5 * (data(numX - 1, -1) + data(numX, 0));
            data(numX, numY) = 0.5 * (data(numX, numY - 1) + data(numX - 1, numY - 1));
            

        }

        void Solver::diffuse(int type, int numX, int numY, Glb::GridData2d& x, Glb::GridData2d& x0, float diff, float dt) {
            float a = dt * diff * numX * numY;
            // ��˹���¶�
            for (int iteration = 100; iteration > 0; iteration--) {
                FOR_EACH_CELL{
                    x(i,j) = (x0(i,j) + a * (x(i + 1,j) + x(i - 1,j) + x(i,j + 1) + x(i,j - 1))) / (1 + 4 * a);
                }
            }
            setbnd(type, numX, numY, x);
        }
        */
    }
}
