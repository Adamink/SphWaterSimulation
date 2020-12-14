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
		// Stores liquid nodes, wall nodes and rigidbody nodes together.
		vector<Particle> nodes;
		Wheel rigidbody;
		vector<Vec3> old_u;
		vector<int> idx_table[l + 1][m + 1][h + 1]; //一张索引表，存放了每个网格内部包含点的索引 \
			                                        //| A index table, which stored all the index of particles in corresponding grids.
		int m_grid_num[l + 1][m + 1][h + 1];
		int m_step; //当前已运行步数 | The step number
		int t_max;//最大的迭代计算步数 | The max steps this programm run.
		int m_ratio_h_dx; // = h/dx
		double m_h, m_dt, m_lambda, m_stiffness, m_KK, m_gamma;
		double dx, dy, dz;
		double m_Kpoly6;
		double a_Kpoly6;
		double vis_Kpoly6;
		double Vmax; //流场中的最大速度 | The max velocity of the fluid region.
		double radius; //绘制图像时的小球半径 | The radius of sphere while drawing the picture.
		double m_mass; //初始流体质量 | The initial mass of fluid particles.
		double rho_L; //初始的流体密度 | The initial density of fluid.
		double rho_R; //初始的刚体密度 | The initial density of rigid body.
		int total_num; //初始粒子总数为0 | The initial number of all the particles is 0.
		double init_pressure; //初始压强 | The initial pressure of particles.
		double vis0; //初始粘度 | The initial viscosity of particles.
		std::string boundary_band; //including 6 char, denotes the boundary of Left/Right/Forward/Back/Down/Up\
								   // P: Periodical; W: Wall.

		int if_dump; //是否要在运行过程中输出文件 | If we need to dump density files while we are running the code.
		int Nwri; // 每隔多少步输出一次文件 | Every "Nwri" steps passed, a density file will be dumped.

		std::vector<Vec3> vertex_localcoord; //存储顶点坐标，坐标系是本地坐标系，即坐标原点在螺旋桨的中心轴上
										//| Store vertex coordinates: A local coordiante, which means the coordinates center is the wheel center
		std::vector<vector<int> > faces; //存储面元信息，使用顶点索引 | Store faces, with vertex index information.
		Vec3 translate_vertex;

		Sph():
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
			vis0(0.000),
			Vmax(0.0),
			if_dump(3),
			Nwri(5),
			m_step(0),
			t_max(1000){
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

			init();
		}

		void init();

		void compute_press_accelerate(int k);

		void compute_vis_accelerate(int k);

		void Compute_Rho(int k);

		void setWallNodes();
		void setRigidBodyMesh();
		void setRigidBodyNodes();
		void setLiquidNodes();

		void addNode(int i, int j, int r, std::string status);

		void dump_file(string file_name);
		void dump_obj_file(string file_name);
		void dump_ply_file(string file_name);
		void dump_wheel_ply(string file_name);
		bool read_file(string file_name);

		void drawing_liquid(void);

		double comp_error();

		void record_velocity();

		void compute(void);

		void step();
		void OnlyReadFileStep();

		int getstep(){ return m_step; }
		int getMaxstep(){ return t_max; }

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