#include "sph.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace SPH
{
	void Sph::add_node(int i, int j, int r, char status){
		if(m_grid_num[i][j][r] == 1){
			if(status != 'L'){
				if(!Nodes[idx_table[i][j][r][0]].isRigidBody() && status == 'R'){
					m_rigidbody.particle_indexes.push_back(idx_table[i][j][r][0]);
				}
				Nodes[idx_table[i][j][r][0]].set_status(status);
			}
			return;
		}
		Particle temp_node(rho_L, p0, g_accelerate);

		idx_table[i][j][r].push_back(total_num++); //Push index into index table
		temp_node.setplace(i * dx, j * dy, r * dz);
		temp_node.set_status(status);
		Nodes.push_back(temp_node);//Push new node into node list.
		if(temp_node.isRigidBody()){ m_rigidbody.particle_indexes.push_back(Nodes.size() - 1); }
		m_grid_num[i][j][r] = 1;
	}

	void Sph::Set_BoundaryNodes(){
		//Special Zone: Wall
		if(boundary_band[0] == 'W'){ //Left
			for(int j = 0; j <= m; j++)
				for(int r = 0; r <= h; r++){
					add_node(0, j, r, 'W'); add_node(1, j, r, 'W');
				}
		}
		if(boundary_band[1] == 'W'){ //Right
			for(int j = 0; j <= m; j++)
				for(int r = 0; r <= h; r++){
					add_node(l, j, r, 'W'); add_node(l - 1, j, r, 'W');
				}
		}
		if(boundary_band[2] == 'W'){ //Forward
			for(int i = 0; i <= l; i++)
				for(int r = 0; r <= h; r++){
					add_node(i, 0, r, 'W'); add_node(i, 1, r, 'W');
				}
		}
		if(boundary_band[3] == 'W'){ //Back
			for(int i = 0; i <= l; i++)
				for(int r = 0; r <= h; r++){
					add_node(i, m, r, 'W'); add_node(i, m - 1, r, 'W');
				}
		}
		if(boundary_band[4] == 'W'){ //Down
			for(int i = 0; i <= l; i++)
				for(int j = 0; j <= m; j++){
					add_node(i, j, 0, 'W'); add_node(i, j, 1, 'W');
				}
		}
		if(boundary_band[5] == 'W'){ //Up
			for(int i = 0; i <= l; i++)
				for(int j = 0; j <= m; j++){
					add_node(i, j, h, 'W'); add_node(i, j, h - 1, 'W');
				}
		}
	}

	void Sph::set_mesh_for_RigidBody(){
		//For wheel, get its mesh:
		double down_height = 2 * dz, up_height = m_rigidbody.m_wheel_height + 2 * dz;
		double r_in = m_rigidbody.m_wheel_radius_insize, r_out = m_rigidbody.m_wheel_radius_outsize;
		vertex_localcoord.push_back(Vec3(0.0f, 0.0f, down_height));//The 2*dy is because the boundary layer.
		vertex_localcoord.push_back(Vec3(0.0f, 0.0f, up_height));//The 2*dy is because the boundary layer.
		double dtheta = 3.1415926 / m_rigidbody.m_leafnum;
		for(int k = 1; k <= m_rigidbody.m_leafnum; k++){
			//step 1: Add 8 new vertex into the list
			int idx_8[8] = { 8 * k - 6, 8 * k - 5, 8 * k - 4, 8 * k - 3, 8 * k - 2, 8 * k - 1, 8 * k, 8 * k + 1 };
			double cos_thetaminor = cos((2 * k - 2) * dtheta), sin_thetaminor = sin((2 * k - 2) * dtheta);
			double cos_thetabigger = cos((2 * k - 1) * dtheta), sin_thetabigger = sin((2 * k - 1) * dtheta);

			vertex_localcoord.push_back(Vec3(r_in * cos_thetaminor, r_in * sin_thetaminor, down_height)); // idx_8[0]
			vertex_localcoord.push_back(Vec3(r_out * cos_thetaminor, r_out * sin_thetaminor, down_height)); // idx_8[1]
			vertex_localcoord.push_back(Vec3(r_out * cos_thetabigger, r_out * sin_thetabigger, down_height)); // idx_8[2]
			vertex_localcoord.push_back(Vec3(r_in * cos_thetabigger, r_in * sin_thetabigger, down_height)); // idx_8[3]

			vertex_localcoord.push_back(Vec3(r_in * cos_thetaminor, r_in * sin_thetaminor, up_height)); // idx_8[4]
			vertex_localcoord.push_back(Vec3(r_out * cos_thetaminor, r_out * sin_thetaminor, up_height)); // idx_8[5]
			vertex_localcoord.push_back(Vec3(r_out * cos_thetabigger, r_out * sin_thetabigger, up_height)); // idx_8[6]
			vertex_localcoord.push_back(Vec3(r_in * cos_thetabigger, r_in * sin_thetabigger, up_height)); // idx_8[7]

			//step 2: Create faces: 4 triangles + 4 Quads = 12 * triangles
			faces.push_back(std::vector<int>{0, idx_8[2], idx_8[1]}); //Trigngle #1
			faces.push_back(std::vector<int>{1, idx_8[5], idx_8[6]}); //Trigngle #2

			faces.push_back(std::vector<int>{idx_8[1], idx_8[2], idx_8[6]}); //Quad #0.5
			faces.push_back(std::vector<int>{idx_8[1], idx_8[6], idx_8[5]}); //Quad #1
			faces.push_back(std::vector<int>{idx_8[2], idx_8[3], idx_8[7]}); //Quad #1.5
			faces.push_back(std::vector<int>{idx_8[2], idx_8[7], idx_8[6]}); //Quad #2
			faces.push_back(std::vector<int>{idx_8[0], idx_8[1], idx_8[5]}); //Quad #2.5
			faces.push_back(std::vector<int>{idx_8[0], idx_8[5], idx_8[4]}); //Quad #3

			int idx_next0 = idx_8[0] + 8, idx_next4 = idx_8[4] + 8;
			if(idx_next4 >= 2 + 8 * m_rigidbody.m_leafnum){
				idx_next0 -= m_rigidbody.m_leafnum * 8;
				idx_next4 -= m_rigidbody.m_leafnum * 8;
			}
			faces.push_back(std::vector<int>{0, idx_next0, idx_8[3]}); //Trigngle #3
			faces.push_back(std::vector<int>{1, idx_8[7], idx_next4}); //Trigngle #4

			faces.push_back(std::vector<int>{idx_8[3], idx_next0, idx_next4}); //Quad #3.5
			faces.push_back(std::vector<int>{idx_8[3], idx_next4, idx_8[7]}); //Quad #4
		}
	}

	void Sph::initialization(){
		Particle temp_node(rho_L, p0, g_accelerate);
		int i = 0, j = 0, r = 0;
		std::cout << "[Init]: Intitialization begin" << std::endl;

		//How to set boundary?
		// if m = 30 (in y direction), then the Boundary is set in idx = 0, 30, and the(1, 29) is fluid.
		// 1. 'W' Points: Wall
		Set_BoundaryNodes();

		// 2. 'R' Points: Rigid body
		set_mesh_for_RigidBody();
		translate_vertex = m_rigidbody.position; translate_vertex.setz(0.0);
		for(i = 0; i <= l; i++){
			for(j = 0; j <= m; j++){
				for(r = 0; r <= h; r++){
					if(r * dz > 2 * dz + m_rigidbody.m_wheel_height){ continue; } //The 2*dy is because the boundary layer.
					if(r * dz < 2 * dz){ continue; }
					double local_x = i * dx - m_rigidbody.position.getx();
					double local_y = j * dy - m_rigidbody.position.gety();
					double r_norm2 = local_x * local_x + local_y * local_y;
					if(r_norm2 > m_rigidbody.m_wheel_radius_outsize * m_rigidbody.m_wheel_radius_outsize){ continue; }
					else if(r_norm2 <= m_rigidbody.m_wheel_radius_insize * m_rigidbody.m_wheel_radius_insize){
						add_node(i, j, r, 'R');
					}
					else{
						double dtheta = 3.1415926 / m_rigidbody.m_leafnum;
						double r_norm = sqrt(r_norm2);
						double costheta = local_x / r_norm, sintheta = local_y / r_norm;
						double theta = acos(costheta);
						if(sintheta < 0.0){ theta = 2 * 3.1415926 - theta; }
						int dnum = (int)((theta + 2 * 3.1415926) / dtheta);
						if(dnum % 2 == 0){
							add_node(i, j, r, 'R');
						}
					}
				}
			}
		}

		// 3. 'L' Points: Liquid
		for(i = 0; i <= l; i++){
			for(j = 0; j < m / 3; j++){
				for(r = 0; r < 3 * h / 4; r++){
					add_node(i, j, r, 'L');
				}
			}
		}

		std::cout << "The total number of particles: " << total_num << std::endl;
		old_u = std::vector<Vec3>(total_num, Vec3(0.0, 0.0, 0.0));
	}

	void Sph::record_velocity(){
		int k = 0;
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			old_u[k] = Nodes[k].velocity;
		}
	}

	double Sph::comp_error(){
		int k = 0;
		double final_error = 0.0;

#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			final_error += (old_u[k] - Nodes[k].velocity).sqrNorm();
		}
		return sqrt(final_error);
	}

	void Sph::compute_vis_accelerate(int k){
		double Vb = 0.0;
		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if(!Nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;

		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if(neighbor_ID == k){ continue; }

			double r_sqnorm = (Nodes[k].position - Nodes[neighbor_ID].position).sqrNorm();
			if(r_sqnorm > m_h * m_h){ continue; }

			Vec3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(Nodes[k].position, Nodes[neighbor_ID].position);
			double u_dot_r = (Nodes[k].velocity - Nodes[neighbor_ID].velocity).dot(Nodes[k].position - Nodes[neighbor_ID].position);
			if(Nodes[neighbor_ID].isFluid()){
				//fluid - fluid
				Nodes[k].force_vis += Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (m_mass * m_mass / Nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
			}
			else{
				//fluid - obstacle
				double mass_rigid = m_mass;
				Vec3 temp_force = Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (mass_rigid * m_mass / Nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
				Nodes[k].force_vis += temp_force;
				Nodes[neighbor_ID].force_vis += -temp_force;
			}
		}
	}

	void Sph::compute_press_accelerate(int k){
		if(Nodes[k].p < 0.0){ return; }

		double Vb = 0.0;
		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if(!Nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;

		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];

			if(neighbor_ID == k){ continue; }
			if((Nodes[k].position - Nodes[neighbor_ID].position).sqrNorm() > m_h * m_h){ continue; }

			Vec3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(Nodes[k].position, Nodes[neighbor_ID].position);
			if(Nodes[neighbor_ID].isWall()){
				//fluid - Obstacle
				double mass_rigid = m_mass;
				Vec3 temp_force = Gradient_Wij * (-2.0 * m_mass * mass_rigid * Nodes[k].p / Nodes[k].rho / Nodes[k].rho);
				Nodes[k].force_press += temp_force;
				Nodes[neighbor_ID].force_press += -temp_force;
			}
			else{
				//fluid - fluid
				Nodes[k].force_press += Gradient_Wij * (-m_mass * m_mass * \
					(Nodes[k].p / Nodes[k].rho / Nodes[k].rho + \
						Nodes[neighbor_ID].p / Nodes[neighbor_ID].rho / Nodes[neighbor_ID].rho));
			}
		}
	}

	void Sph::Compute_Rho(int k){
		if(Nodes[k].isWall()){ return; }
		Nodes[k].rho = 0.0;
		double Vb = 0.0;
		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if(!Nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;
		for(int i = 0; i < Nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if((Nodes[k].position - Nodes[neighbor_ID].position).sqrNorm() > m_h * m_h){ continue; }

			Vec3 Grad_W = Gradient_Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
			if(!Nodes[neighbor_ID].isFluid()){
				double mass_rigid = m_mass;
				Nodes[k].rho += mass_rigid * Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
				Nodes[k].rho += mass_rigid * m_dt * (Nodes[k].velocity - Nodes[neighbor_ID].velocity).dot(Grad_W);
			}
			else{
				Nodes[k].rho += m_mass * Muller03Kernel_Basic(Nodes[k].position, Nodes[neighbor_ID].position);
				Nodes[k].rho += m_mass * m_dt * (Nodes[k].velocity - Nodes[neighbor_ID].velocity).dot(Grad_W);
			}
		}
	}

	void Sph::compute(void){
		int i = 0, j = 0, r = 0, k = 0;

		//step0:更新表：
#pragma omp parallel for private(i,j,r)
		for(i = 0; i <= l; i++){
			for(j = 0; j <= m; j++){
				for(r = 0; r <= h; r++){
					idx_table[i][j][r].clear();
					m_grid_num[i][j][r] = 0;
				}
			}
		}
		for(k = 0; k < total_num; k++){
			i = (int)(Nodes[k].position.getx() / dx);
			j = (int)(Nodes[k].position.gety() / dy);
			r = (int)(Nodes[k].position.getz() / dz);
			if(i < 0 || i > l){ std::cerr << "A particle flied out from x direction:" << Nodes[k].status << std::endl; continue; }
			if(j < 0 || j > m){ std::cerr << "A particle flied out from y direction:" << Nodes[k].status << std::endl; continue; }
			if(r < 0 || r > h){ std::cerr << "A particle flied out from z direction:" << Nodes[k].status << std::endl; continue; }
			idx_table[i][j][r].push_back(k);
			m_grid_num[i][j][r] += 1;
		}
		for(k = 0; k < total_num; k++){
			Nodes[k].force_vis = Vec3(0.0, 0.0, 0.0);
			Nodes[k].neighbor_index.clear();
			int comp_i = (int)(Nodes[k].position.getx() / dx);
			int comp_j = (int)(Nodes[k].position.gety() / dy);
			int comp_r = (int)(Nodes[k].position.getz() / dz);
			for(i = ((comp_i >= m_ratio_h_dx) ? (comp_i - m_ratio_h_dx) : 0); i <= ((comp_i <= l - m_ratio_h_dx) ? (comp_i + m_ratio_h_dx) : l); i++){
				for(j = ((comp_j >= m_ratio_h_dx) ? (comp_j - m_ratio_h_dx) : 0); j <= ((comp_j <= m - m_ratio_h_dx) ? (comp_j + m_ratio_h_dx) : m); j++){
					for(r = ((comp_r >= m_ratio_h_dx) ? (comp_r - m_ratio_h_dx) : 0); r <= ((comp_r <= h - m_ratio_h_dx) ? (comp_r + m_ratio_h_dx) : h); r++){
						for(int comp_index = 0; comp_index < m_grid_num[i][j][r]; comp_index++){
							int neighbor_ID = idx_table[i][j][r][comp_index];
							double r_sqnorm = (Nodes[k].position - Nodes[neighbor_ID].position).sqrNorm();
							if(r_sqnorm <= m_h * m_h){ Nodes[k].neighbor_index.push_back(neighbor_ID); }
						}
					}
				}
			}
		}

		//step1:计算粒子的Non-Pressure加速度：
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			//Nodes[k].force_vis = Vec3(0.0, 0.0, 0.0); execute in step 0.
			//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
			if(Nodes[k].isWall()){ continue; }
			compute_vis_accelerate(k);
		}

		//step2:更新粒子的速度
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			if(Nodes[k].isWall()){ continue; }
			Nodes[k].velocity += (Nodes[k].force_vis / m_mass + Nodes[k].acc_ext) * m_dt;
		}

		int iter_num = 0;
		double rho_error = 0.0;
		while(++iter_num <= 10){
			//SteVec3: 根据光滑核函数计算粒子的插值密度：
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				Compute_Rho(k);
			}

			//Step4: 计算粒子的压强：
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				Nodes[k].force_press = Vec3(0.0, 0.0, 0.0);
				Nodes[k].p = m_KK * (pow(Nodes[k].rho / rho_L, (int)m_gamma) - 1.0);
				if(Nodes[k].p < 0.0){ Nodes[k].p = 0.0; }
			}

			//step5:计算粒子的Pressure加速度：
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				//Nodes[k].acc_press = Vec3(0.0, 0.0, 0.0); execute in step 4.
				//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
				if(!Nodes[k].isFluid()){ continue; }
				compute_press_accelerate(k);
			}

			//step6:更新粒子的速度
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				if(!Nodes[k].isFluid()){ continue; }
				Nodes[k].velocity += Nodes[k].force_press * (m_dt / m_mass);
			}

			//Check if we can break out the iteration:
			rho_error = 0.0;
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				rho_error += (Nodes[k].rho - rho_L) * (Nodes[k].rho - rho_L);
			}
			rho_error = sqrt(rho_error);
			if(rho_error < 0.1){ break; }
		}
		std::cerr << "rho_error: " << rho_error << std::endl;

		//step6.5:更新刚体运动
		m_rigidbody.update(m_dt);
		double dTheta_dt = m_rigidbody.swirl_velocity;
		double dTheta = dTheta_dt * m_dt;
		double cosdTheta = cos(dTheta), sindTheta = sin(dTheta);
#pragma omp parallel for private(k,i)
		for(i = 0; i < m_rigidbody.particle_indexes.size(); i++){
			k = m_rigidbody.particle_indexes[i];
			double local_x = Nodes[k].position.getx() - m_rigidbody.position.getx();
			double local_y = Nodes[k].position.gety() - m_rigidbody.position.gety();
			double r_norm = sqrt(local_x * local_x + local_y * local_y);
			double costheta = local_x / r_norm, sintheta = local_y / r_norm;
			//update:
			Nodes[k].setplace((costheta * cosdTheta - sintheta * sindTheta) * r_norm + m_rigidbody.position.getx(),
				(sintheta * cosdTheta + costheta * sindTheta) * r_norm + m_rigidbody.position.gety(),
				Nodes[k].position.getz());
		}
		Mat3 RotationMat(std::vector<Vec3>{Vec3(cosdTheta, -sindTheta, 0.0),
			Vec3(sindTheta, cosdTheta, 0.0),
			Vec3(0.0, 0.0, 1.0)});
		for(i = 0; i < vertex_localcoord.size(); i++){
			vertex_localcoord[i] = RotationMat.multiply(vertex_localcoord[i]);
		}


		//step7:更新粒子的位置：
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			if(!Nodes[k].isFluid()){ continue; }
			Nodes[k].position += Nodes[k].velocity * m_dt;
			if(boundary_band[0] == 'P'){ //Left
				if(Nodes[k].position.getx() < 0){ Nodes[k].position += Vec3(Lx, 0.0, 0.0); }
			}
			if(boundary_band[1] == 'P'){ //Right
				if(Nodes[k].position.getx() > Lx){ Nodes[k].position += Vec3(-Lx, 0.0, 0.0); }
			}
			if(boundary_band[2] == 'P'){ //Forward
				if(Nodes[k].position.gety() < 0){ Nodes[k].position += Vec3(0.0, Ly, 0.0); }
			}
			if(boundary_band[3] == 'P'){ //Back
				if(Nodes[k].position.gety() > Ly){ Nodes[k].position += Vec3(0.0, -Ly, 0.0); }
			}
			if(boundary_band[4] == 'P'){ //Down
				if(Nodes[k].position.getz() < 0){ Nodes[k].position += Vec3(0.0, 0.0, Lz); }
			}
			if(boundary_band[5] == 'P'){ //Up
				if(Nodes[k].position.getz() > Lz){ Nodes[k].position += Vec3(0.0, 0.0, -Lz); }
			}
		}

		Vmax = 0.0;
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			double u_maxAbsCoord = Nodes[k].velocity.maxAbsCoord();
			Vmax = Vmax > u_maxAbsCoord ? Vmax : u_maxAbsCoord;
		}
		//update m_dt:
		m_dt = m_lambda * m_h / max(Vmax, sqrt(m_stiffness));
		m_dt = m_lambda * m_h / sqrt(m_stiffness);
		std::cout << std::endl;

		return;
	}

	void Sph::dump_file(string file_name){
		//本函数用于输出可以被软件OVITO可视化的文件 | This function can be used to dump a file which can be visualized by software OVITO
		ofstream fout(file_name);
		if(!fout){
			std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
			exit(-1);
		}

		fout << "Number of particles = " << total_num << std::endl;
		fout << "A = 1 Angstrom (basic length-scale)" << std::endl;

		fout << "H0(1,1) = " << Lx * 100 << " A" << std::endl;
		fout << "H0(1,2) = " << 0 << " A" << std::endl;
		fout << "H0(1,3) = " << 0 << " A" << std::endl;

		fout << "H0(2,1) = " << 0 << " A" << std::endl;
		fout << "H0(2,2) = " << Ly * 100 << " A" << std::endl;
		fout << "H0(2,3) = " << 0 << " A" << std::endl;

		fout << "H0(3,1) = " << 0 << " A" << std::endl;
		fout << "H0(3,2) = " << 0 << " A" << std::endl;
		fout << "H0(3,3) = " << Lz * 100 << " A" << std::endl;

		fout << ".NO_VELOCITY." << std::endl;
		fout << "entry_count = " << 3 << std::endl;

		for(int k = 0; k < total_num; k++){
			fout << m_mass << std::endl; // mass
			fout << Nodes[k].status << std::endl; // element type
			fout << Nodes[k].position.getx() * 100 << " " << Nodes[k].position.gety() * 100 << " " << Nodes[k].position.getz() * 100 << std::endl; // Info of position.
		}

		fout.close();

		return;
	}

	void Sph::dump_wheel_ply(string file_name){
		//本函数用于输出ply文件以显示刚体的网格信息 | This function can be used to dump a ply file for rigidbody's mesh
		ofstream fout(file_name);
		if(!fout){
			std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
			exit(-1);
		}
		fout << "ply" << std::endl;
		fout << "format ascii 1.0" << std::endl;
		fout << "comment made by anonymous" << std::endl;
		fout << "comment this file is a wheel" << std::endl;
		fout << "element vertex " << vertex_localcoord.size() << std::endl;
		fout << "property float32 x" << std::endl;
		fout << "property float32 y" << std::endl;
		fout << "property float32 z" << std::endl;
		fout << "element face " << faces.size() << std::endl;
		fout << "property list uint8 int32 vertex_index" << std::endl;
		fout << "end_header" << std::endl;
		//输出顶点坐标 | dump vertex coordinate into file
		for(int k = 0; k < vertex_localcoord.size(); k++){
			Vec3 single_vertex_worldcoord = vertex_localcoord[k] + translate_vertex;
			fout << single_vertex_worldcoord.getx() << " " << single_vertex_worldcoord.gety() << " " << \
				single_vertex_worldcoord.getz() << std::endl;
		}
		//输出面元的顶点索引 | dump the vertex index of faces into file
		for(int k = 0; k < faces.size(); k++){
			fout << "3 ";
			for(int rr = 0; rr < faces[k].size(); rr++){
				fout << faces[k][rr] << " ";
			}
			fout << std::endl;
		}
	}


	void Sph::dump_ply_file(string file_name){
		stringstream ss;
		ofstream fout(file_name);
		if(!fout){
			std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
			exit(-1);
		}
		int cnt = 0;
		for(int k = 0; k < total_num; k++){
			if(Nodes[k].status == 'W' || Nodes[k].status == 'R') continue; // note Wall or Rigid Body
			cnt++;
			ss << Nodes[k].position.getx() << " " << Nodes[k].position.gety() << " " << Nodes[k].position.getz() << endl;
		}
		fout << "ply" << endl;
		fout << "format ascii 1.0" << endl;
		fout << "element vertex " << cnt << endl;
		fout << "property float x" << endl;
		fout << "property float y" << endl;
		fout << "property float z" << endl;
		fout << "end_header" << endl;
		fout << ss.str();
	}

	bool Sph::read_file(string file_name){
		//本函数用于读取已经之前本程序输出过的文件
		ifstream fin(file_name);
		string get_str;
		if(!fin){
			std::cout << "(in_put function)Error! Can not read this file." << std::endl;
			return false;
		}
		for(int k = 0; k < 13; k++){
			//line1: Number of particles = {$total_num}
			//line2: A = 1 Angstrom (basic length-scale)
			//line3: H0(1,1) = {$Lx * 100} A
			//line4: H0(1,2) = 0 A
			//line5: H0(1,3) = 0 A
			//line6: H0(2,1) = 0 A
			//line7: H0(2,2) = {$Ly * 100} A
			//line8: H0(2,3) = 0 A
			//line9: H0(3,1) = 0 A
			//line10: H0(3,2) = 0 A
			//line11: H0(3,3) = {$L * 100} A
			//line12: .NO_VELOCITY.
			//line13: entry_count = 3
			getline(fin, get_str);
		}
		int k = 0;
		double in_x, in_y, in_z;
		char status;
		while(k != total_num){
			fin >> m_mass;
			fin >> status;
			Nodes[k].set_status(status);
			fin >> in_x >> in_y >> in_z;
			Nodes[k].setplace(in_x / 100.0, in_y / 100.0, in_z / 100.0);
			k++;
		}
		std::cout << "m_step: " << m_step << ", total_num: " << total_num << std::endl;
		return true;
	}

	void Sph::drawing_liquid(void){
		GLint i = 0, j = 0, r = 0;

		glColor3f(0.0, 1.0, 0.0);
		//绘制流场外围
		glBegin(GL_LINE_STRIP);
		glVertex3i(0, 0, 0);
		glVertex3i(l, 0, 0);
		glVertex3i(l, m, 0);
		glVertex3i(0, m, 0);
		glVertex3i(0, 0, 0);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex3i(0, 0, h);
		glVertex3i(l, 0, h);
		glVertex3i(l, m, h);
		glVertex3i(0, m, h);
		glVertex3i(0, 0, h);
		glEnd();
		glBegin(GL_LINES);
		glVertex3i(0, 0, 0);
		glVertex3i(0, 0, h);
		glVertex3i(l, 0, 0);
		glVertex3i(l, 0, h);
		glVertex3i(0, m, 0);
		glVertex3i(0, m, h);
		glVertex3i(l, m, 0);
		glVertex3i(l, m, h);
		glEnd();
		//绘制流场外围结束

		GLfloat* vertices = nullptr;
		vertices = (GLfloat*)malloc(3 * total_num * (sizeof(GLfloat)));
		for(i = total_num - 1; i >= 0; i--){
			glColor3f(0.5f, 0.5f, 1.0f);
			if(Nodes[i].isWall()){
				continue;
			}
			if(Nodes[i].isRigidBody()){
				glColor3f(1.0f, 1.0f, 0.0f);
			}
			//if (Nodes[i].position[0] < 0 || Nodes[i].position[0] > Lx) { continue; }
			//if (Nodes[i].position[1] < 0 || Nodes[i].position[1] > Ly) { continue; }
			//if (Nodes[i].position[2] < 0 || Nodes[i].position[2] > Lz) { continue; }
			glPushMatrix();
			glTranslated(Nodes[i].position.getx() / dx, Nodes[i].position.gety() / dy, Nodes[i].position.getz() / dz);
			glutSolidSphere(radius, 6, 6);
			glPopMatrix();

			//glVertex3f(Nodes[i].position[0] / dx, Nodes[i].position[1] / dy, Nodes[i].position[2] / dz);

			//*(vertices + 3 * i + 0) = Nodes[i].position[0] / dx;
			//*(vertices + 3 * i + 1) = Nodes[i].position[1] / dy;
			//*(vertices + 3 * i + 2) = Nodes[i].position[2] / dz;
		}

		//if(m_openglcontext != nullptr)
			//m_openglcontext->VBO_bind<GLfloat>(vertices);

		return;
	}

	void Sph::step(){
		double final_error = 1.0;
		record_velocity();
		drawing_liquid();
		compute();
		cout << "step " << m_step << endl;
		if(m_step++ % Nwri == 0){
			final_error = comp_error();
			std::cout << "=====================" << std::endl;
			std::cout << "Step: " << m_step << ";  error: " << final_error << std::endl;
			// char str_step[10];
			// itoa(m_step, str_step, 10);
			stringstream ss;
			ss << std::setfill('0') << std::setw(4) << m_step / Nwri;
			if(if_dump == 1){ dump_file(string("data") + ss.str() + ".cfg"); }
			else if(if_dump == 3){ dump_ply_file("../output/liquid/" + ss.str() + ".ply"); }
			// in all cases dump wheel 
			dump_wheel_ply("../output/wheel/" + ss.str() + ".ply");

		}
	}

	void Sph::OnlyReadFileStep(){
		drawing_liquid();
		if(m_step++ % Nwri == 0){
			// char str_step[10];
			// itoa(m_step, str_step, 10);
			stringstream ss;
			ss << m_step;
			if(!read_file(string("data") + ss.str() + ".cfg")){
				m_step = 0;
			}
		}
	}
}