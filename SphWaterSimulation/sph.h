#ifndef SPH_H_
#define SPH_H_

#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include <GL/glut.h>

#include "particle.h"
#include "vec3.h"
#include "rigid_body.h"
#include "constants.h"

using namespace std;

namespace SPH
{
	// Define the size of bounding box of simulation with x in [0, x_bound], y in [0, y_bound], z in [0, z_bound].
	static const double x_bound = sph_const::kXBound, y_bound = sph_const::kYBound, z_bound = sph_const::kZBound;

	// Define the number of particles, and in total l x m x h particles are inited on 3d-grids.
	static const int l = sph_const::kXParticleNum, m = sph_const::kYParticleNum, h = sph_const::kZParticleNum;
	static const Vec3 g_gravity = Vec3(0.0, 0.0, -sph_const::kGravity);

	class Sph{
	public:
		bool use_rigid_body;
		bool if_visualize;
		bool read_from_file;
		// Store liquid nodes, wall nodes and rigidbody nodes together.
		vector<Particle> nodes;
		// Store all the index of particles in corresponding grids
		vector<int> idx_table[l + 1][m + 1][h + 1];
		Wheel rigidbody;
		vector<Vec3> prevVelocity;
		int m_grid_num[l + 1][m + 1][h + 1];
		int m_step;  // Current step number
		int m_ratio_h_dx; // = h/dx
		double m_h, m_dt, m_lambda, m_stiffness, m_KK, m_gamma;
		double dx, dy, dz;
		double m_Kpoly6;
		double a_Kpoly6;
		double vis_Kpoly6;
		double Vmax; //流场中的最大速度 | The max velocity of the fluid region.
		double radius; //绘制图像时的小球半径 | The radius of sphere while draw the picture.
		double m_mass; //初始流体质量 | The initial mass of fluid particles.
		double rho_L; //初始的流体密度 | The initial density of fluid.
		double rho_R; //初始的刚体密度 | The initial density of rigid body.
		int total_num; //初始粒子总数为0 | The initial number of all the particles is 0.
		double init_pressure; //初始压强 | The initial pressure of particles.
		double viscosity; //初始粘度 | The initial viscosity of particles.

		// Including 6 char, each denoting the boundary type of Left/Right/Forward/Back/Down/Up,
		// where 'P' stands for Periodical and 'W' stands for Wall
		std::string boundary_band;

		int dump_file_interval; // 每隔多少步输出一次文件 | Every "dump_file_interval" steps passed, a density file will be dumped.

		std::vector<Vec3> vertex_localcoord; //存储顶点坐标，坐标系是本地坐标系，即坐标原点在螺旋桨的中心轴上
										//| Store vertex coordinates: A local coordiante, which means the coordinates center is the wheel center
		std::vector<vector<int> > faces; //存储面元信息，使用顶点索引 | Store faces, with vertex index information.
		Vec3 translate_vertex;

		Sph():
			use_rigid_body(program_const::kUseRigidBody),
			if_visualize(program_const::IF_VISUALIZE),
			read_from_file(program_const::READ_FROM_FILES),
			m_ratio_h_dx(2),
			dx(x_bound / l),
			dy(y_bound / m),
			dz(z_bound / h),
			m_lambda(0.4),
			m_stiffness(1000.0),
			m_mass(1.0),
			m_gamma(1.0),
			radius(0.4),
			total_num(0),
			init_pressure(0.0),
			viscosity(sph_const::kViscosity),
			Vmax(0.0),
			dump_file_interval(program_const::kDumpFileInterval),
			m_step(0){
			m_h = m_ratio_h_dx * dx;
			rho_L = m_mass / dx / dy / dz;
			rho_R = 300.0;
			m_dt = m_lambda * m_h / sqrt(m_stiffness);
			m_KK = m_stiffness * rho_L / m_gamma;
			m_Kpoly6 = 315 / ((64 * math_const::PI) * pow(m_h, 9));
			a_Kpoly6 = 45 / (math_const::PI * pow(m_h, 6));
			vis_Kpoly6 = 15.0 / (2.0 * math_const::PI * pow(m_h, 3));
			//rigidbody = Sphere(rho_R, y_bound / 6.0, Vec3(x_bound / 2.0, 2.0 * y_bound / 3.0, y_bound / 6.0 + 2 * dx));
			rigidbody = Wheel(rho_R, x_bound / 6.0, x_bound / 12.0, z_bound / 5.0, Vec3(x_bound / 2.0, y_bound / 2.0, 0.0), 4);
			for(int i = 0; i <= l; i++){
				for(int j = 0; j <= m; j++){
					for(int r = 0; r <= h; r++){
						m_grid_num[i][j][r] = 0;
					}
				}
			}

			initNodes();
		}

		void initNodes();

		void computePressureAccelerate(int k);

		void computeVisAccerlerate(int k);

		void computeDensity(int k);

		void setWallNodes();
		void setRigidBodyMesh();
		void setRigidBodyNodes();
		void setLiquidNodes();

		void updateRigidBody();
		void addNode(int i, int j, int r, std::string status);

		std::string getFilePath(std::string command);
		void dumpLiquidAsCfg(std::string file_path);
		void dumpLiquidAsPly(std::string file_path);
		void dumpRigidBody(std::string file_path);
		bool readLiquidFromCfg(std::string file_path);

		void draw();
		void drawBoundingBox();
		void drawRigidBody();
		void drawLiquid();

		double getError();

		void recordVelocity();

		void compute();

		void step();
		void readFileAndShow();

		int getStep(){ return m_step; }

		double Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point){
			double r_square_norm = (detect_point - origin_point).sqrNorm();
			if(m_h * m_h > r_square_norm){ return m_Kpoly6 * pow(m_h * m_h - r_square_norm, 3); }
			else{ return 0.0; }
		}

		Vec3 Gradient_Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point){
			double r_square_norm = (detect_point - origin_point).sqrNorm();
			if(m_h * m_h > r_square_norm){
				return (detect_point - origin_point) * (6.0 * m_Kpoly6 * pow(m_h * m_h - r_square_norm, 2));
			}
			else{ return Vec3(0.0, 0.0, 0.0); }
		}

		double Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point){
			double r_norm = (origin_point - detect_point).norm();
			if(m_h > r_norm){ return a_Kpoly6 / 3.0 * pow((m_h - r_norm), 3); }
			else{ return 0.0; }
		}

		Vec3 Gradient_Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point){
			double r_norm = (origin_point - detect_point).norm();
			if(r_norm > m_h){ return Vec3(0.0, 0.0, 0.0); }
			return (detect_point - origin_point).normalized() * (a_Kpoly6 * pow(m_h - r_norm, 2));
		}

		Vec3 Gradient_Muller03Kernel_Vis(Vec3 origin_point, Vec3 detect_point){
			double r_norm = (origin_point - detect_point).norm();
			if(r_norm > m_h){ return Vec3(0.0, 0.0, 0.0); }
			return (detect_point - origin_point) * (vis_Kpoly6 * \
				(-1.5 * r_norm / pow(m_h, 3) + 2.0 / m_h / m_h - 0.5 * m_h / pow(r_norm, 3)));
		}
	};
}

#endif