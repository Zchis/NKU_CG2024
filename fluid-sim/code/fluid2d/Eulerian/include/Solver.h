#pragma once
#ifndef __EULERIAN_2D_SOLVER_H__
#define __EULERIAN_2D_SOLVER_H__

#include "Eulerian/include/MACGrid2d.h"
#include "Global.h"

/*
namespace FluidSimulation {
	namespace Eulerian2d {
		class Solver {
		public:
			Solver(MACGrid2d& grid);

			void solve();

		protected:
			MACGrid2d& mGrid;
		};
	}
}

#endif // !__EULER_SOLVER_H__
*/

namespace FluidSimulation {
    namespace Eulerian2d {
        class Solver {
        public:
            Solver(MACGrid2d& grid);

            void solve();

        protected:

            void vel_step(float dt);
            void dens_step(float dt);

            // 1. advection
            // 2. compute external forces
            // 3. projection
            void advect(float dt);

            void computeforces(float dt);

            void project(float dt);

            MACGrid2d& mGrid;
        };
    }
}

#endif // !__EULERIAN_2D_SOLVER_H__