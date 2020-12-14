#include "sph.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "constants.h"

namespace SPH
{
	void Sph::init(){
		Particle temp_node(rho_L, init_pressure, g_gravity);
		// 1. Set "Wall" nodes as the bounding box.
		setWallNodes();
		// 2. Set "RigidBody" nodes as the wheel, as well as the meshes for rendering.
		setRigidBodyMesh();
		setRigidBodyNodes();
		// 3. Set "Liquid" nodes as the water.
		setLiquidNodes();

		std::cout << "The total number of particles: " << total_num << std::endl;
		old_u = std::vector<Vec3>(total_num, Vec3(0.0, 0.0, 0.0));
	}

	// Register a node of different kinds to 'nodes' as well as 'idx_table' during init().
	// If the node is rigidbody
	void Sph::addNode(int i, int j, int r, std::string status){
		// Rigidbody and wall nodes takes the first priority.
		if(m_grid_num[i][j][r] == 1){
			if(status != "Liquid"){
				if(!nodes[idx_table[i][j][r][0]].isRigidBody() && status == "RigidBody"){
					rigidbody.particle_indexes.push_back(idx_table[i][j][r][0]);
				}
				nodes[idx_table[i][j][r][0]].setStatus(status);
			}
			return;
		}
		Particle temp_node(rho_L, init_pressure, g_gravity);
		// Register the index of particle to 'idx_table'.
		idx_table[i][j][r].push_back(total_num++);
		temp_node.setPosition(i * dx, j * dy, r * dz);
		temp_node.setStatus(status);
		nodes.push_back(temp_node);//Push new node into node list.
		if(temp_node.isRigidBody()){ rigidbody.particle_indexes.push_back(nodes.size() - 1); }
		m_grid_num[i][j][r] = 1;
	}

	// Set "Wall" nodes as the bounding box
	// Example: if m = 30 in y direction, then particles with index 0, 30 will be wall nodes,
	// and particles with index from 1 to 29 (inclusive) will be fluid nodes
	void Sph::setWallNodes(){
		for(int j = 0; j <= m; j++)
			for(int r = 0; r <= h; r++){
				addNode(0, j, r, "Wall"); addNode(1, j, r, "Wall");
			}
		for(int j = 0; j <= m; j++)
			for(int r = 0; r <= h; r++){
				addNode(l, j, r, "Wall"); addNode(l - 1, j, r, "Wall");
			}
		for(int i = 0; i <= l; i++)
			for(int r = 0; r <= h; r++){
				addNode(i, 0, r, "Wall"); addNode(i, 1, r, "Wall");
			}
		for(int i = 0; i <= l; i++)
			for(int r = 0; r <= h; r++){
				addNode(i, m, r, "Wall"); addNode(i, m - 1, r, "Wall");
			}
		for(int i = 0; i <= l; i++)
			for(int j = 0; j <= m; j++){
				addNode(i, j, 0, "Wall"); addNode(i, j, 1, "Wall");
			}
		for(int i = 0; i <= l; i++)
			for(int j = 0; j <= m; j++){
				addNode(i, j, h, "Wall"); addNode(i, j, h - 1, "Wall");
			}
	}

	void Sph::setRigidBodyNodes(){
		translate_vertex = rigidbody.position; translate_vertex.setz(0.0);
		int i = 0, j = 0, r = 0;
		for(i = 0; i <= l; i++){
			for(j = 0; j <= m; j++){
				for(r = 0; r <= h; r++){
					if(r * dz > 2 * dz + rigidbody.m_wheel_height){ continue; } //The 2*dy is because the boundary layer.
					if(r * dz < 2 * dz){ continue; }
					double local_x = i * dx - rigidbody.position.getx();
					double local_y = j * dy - rigidbody.position.gety();
					double r_norm2 = local_x * local_x + local_y * local_y;
					if(r_norm2 > rigidbody.m_wheel_radius_outsize * rigidbody.m_wheel_radius_outsize){ continue; }
					else if(r_norm2 <= rigidbody.m_wheel_radius_insize * rigidbody.m_wheel_radius_insize){
						addNode(i, j, r, "RigidBody");
					}
					else{
						double dtheta = math_const::PI / rigidbody.m_leafnum;
						double r_norm = sqrt(r_norm2);
						double costheta = local_x / r_norm, sintheta = local_y / r_norm;
						double theta = acos(costheta);
						if(sintheta < 0.0){ theta = 2 * math_const::PI - theta; }
						int dnum = (int)((theta + 2 * math_const::PI) / dtheta);
						if(dnum % 2 == 0){
							addNode(i, j, r, "RigidBody");
						}
					}
				}
			}
		}
	}


	// 
	void Sph::setRigidBodyMesh(){
		double down_height = 2 * dz, up_height = rigidbody.m_wheel_height + 2 * dz;
		double r_in = rigidbody.m_wheel_radius_insize, r_out = rigidbody.m_wheel_radius_outsize;
		vertex_localcoord.push_back(Vec3(0.0f, 0.0f, down_height));//The 2*dy is because the boundary layer.
		vertex_localcoord.push_back(Vec3(0.0f, 0.0f, up_height));//The 2*dy is because the boundary layer.
		double dtheta = math_const::PI / rigidbody.m_leafnum;
		for(int k = 1; k <= rigidbody.m_leafnum; k++){
			// Step 1: Add 8 new vertex into the list
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

			// Step 2: Create faces: 4 triangles + 4 Quads = 12 * triangles
			faces.push_back(std::vector<int>{0, idx_8[2], idx_8[1]}); //Trigngle #1
			faces.push_back(std::vector<int>{1, idx_8[5], idx_8[6]}); //Trigngle #2

			faces.push_back(std::vector<int>{idx_8[1], idx_8[2], idx_8[6]}); //Quad #0.5
			faces.push_back(std::vector<int>{idx_8[1], idx_8[6], idx_8[5]}); //Quad #1
			faces.push_back(std::vector<int>{idx_8[2], idx_8[3], idx_8[7]}); //Quad #1.5
			faces.push_back(std::vector<int>{idx_8[2], idx_8[7], idx_8[6]}); //Quad #2
			faces.push_back(std::vector<int>{idx_8[0], idx_8[1], idx_8[5]}); //Quad #2.5
			faces.push_back(std::vector<int>{idx_8[0], idx_8[5], idx_8[4]}); //Quad #3

			int idx_next0 = idx_8[0] + 8, idx_next4 = idx_8[4] + 8;
			if(idx_next4 >= 2 + 8 * rigidbody.m_leafnum){
				idx_next0 -= rigidbody.m_leafnum * 8;
				idx_next4 -= rigidbody.m_leafnum * 8;
			}
			faces.push_back(std::vector<int>{0, idx_next0, idx_8[3]}); //Trigngle #3
			faces.push_back(std::vector<int>{1, idx_8[7], idx_next4}); //Trigngle #4

			faces.push_back(std::vector<int>{idx_8[3], idx_next0, idx_next4}); //Quad #3.5
			faces.push_back(std::vector<int>{idx_8[3], idx_next4, idx_8[7]}); //Quad #4
		}
	}

	void Sph::setLiquidNodes(){
		for(int i = 0; i <= l; i++){
			for(int j = 0; j < m / 3; j++){
				for(int r = 0; r < 3 * h / 4; r++){
					addNode(i, j, r, "Liquid");
				}
			}
		}
	}

	void Sph::record_velocity(){
		int k = 0;
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			old_u[k] = nodes[k].velocity;
		}
	}

	double Sph::comp_error(){
		int k = 0;
		double final_error = 0.0;

#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			final_error += (old_u[k] - nodes[k].velocity).sqrNorm();
		}
		return sqrt(final_error);
	}

	void Sph::compute_vis_accelerate(int k){
		double Vb = 0.0;
		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];
			if(!nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;

		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];
			if(neighbor_ID == k){ continue; }

			double r_sqnorm = (nodes[k].position - nodes[neighbor_ID].position).sqrNorm();
			if(r_sqnorm > m_h * m_h){ continue; }

			Vec3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(nodes[k].position, nodes[neighbor_ID].position);
			double u_dot_r = (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(nodes[k].position - nodes[neighbor_ID].position);
			if(nodes[neighbor_ID].isFluid()){
				//fluid - fluid
				nodes[k].force_vis += Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (m_mass * m_mass / nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
			}
			else{
				//fluid - obstacle
				double mass_rigid = m_mass;
				Vec3 temp_force = Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (mass_rigid * m_mass / nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
				nodes[k].force_vis += temp_force;
				nodes[neighbor_ID].force_vis += -temp_force;
			}
		}
	}

	void Sph::compute_press_accelerate(int k){
		if(nodes[k].pressure < 0.0){ return; }

		double Vb = 0.0;
		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];
			if(!nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;

		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];

			if(neighbor_ID == k){ continue; }
			if((nodes[k].position - nodes[neighbor_ID].position).sqrNorm() > m_h * m_h){ continue; }

			Vec3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(nodes[k].position, nodes[neighbor_ID].position);
			if(nodes[neighbor_ID].isWall()){
				//fluid - Obstacle
				double mass_rigid = m_mass;
				Vec3 temp_force = Gradient_Wij * (-2.0 * m_mass * mass_rigid * nodes[k].pressure / nodes[k].rho / nodes[k].rho);
				nodes[k].force_press += temp_force;
				nodes[neighbor_ID].force_press += -temp_force;
			}
			else{
				//fluid - fluid
				nodes[k].force_press += Gradient_Wij * (-m_mass * m_mass * \
					(nodes[k].pressure / nodes[k].rho / nodes[k].rho + \
						nodes[neighbor_ID].pressure / nodes[neighbor_ID].rho / nodes[neighbor_ID].rho));
			}
		}
	}

	void Sph::Compute_Rho(int k){
		if(nodes[k].isWall()){ return; }
		nodes[k].rho = 0.0;
		double Vb = 0.0;
		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];
			if(!nodes[neighbor_ID].isFluid()){ //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
			}
		}
		Vb = 1.0 / Vb;
		for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
			int neighbor_ID = nodes[k].neighbor_index[i];
			if((nodes[k].position - nodes[neighbor_ID].position).sqrNorm() > m_h * m_h){ continue; }

			Vec3 Grad_W = Gradient_Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
			if(!nodes[neighbor_ID].isFluid()){
				double mass_rigid = m_mass;
				nodes[k].rho += mass_rigid * Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
				nodes[k].rho += mass_rigid * m_dt * (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(Grad_W);
			}
			else{
				nodes[k].rho += m_mass * Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
				nodes[k].rho += m_mass * m_dt * (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(Grad_W);
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
			i = (int)(nodes[k].position.getx() / dx);
			j = (int)(nodes[k].position.gety() / dy);
			r = (int)(nodes[k].position.getz() / dz);
			if(i < 0 || i > l){ std::cerr << "A particle flied out from x direction:" << nodes[k].status << std::endl; continue; }
			if(j < 0 || j > m){ std::cerr << "A particle flied out from y direction:" << nodes[k].status << std::endl; continue; }
			if(r < 0 || r > h){ std::cerr << "A particle flied out from z direction:" << nodes[k].status << std::endl; continue; }
			idx_table[i][j][r].push_back(k);
			m_grid_num[i][j][r] += 1;
		}
		for(k = 0; k < total_num; k++){
			nodes[k].force_vis = Vec3(0.0, 0.0, 0.0);
			nodes[k].neighbor_index.clear();
			int comp_i = (int)(nodes[k].position.getx() / dx);
			int comp_j = (int)(nodes[k].position.gety() / dy);
			int comp_r = (int)(nodes[k].position.getz() / dz);
			for(i = ((comp_i >= m_ratio_h_dx) ? (comp_i - m_ratio_h_dx) : 0); i <= ((comp_i <= l - m_ratio_h_dx) ? (comp_i + m_ratio_h_dx) : l); i++){
				for(j = ((comp_j >= m_ratio_h_dx) ? (comp_j - m_ratio_h_dx) : 0); j <= ((comp_j <= m - m_ratio_h_dx) ? (comp_j + m_ratio_h_dx) : m); j++){
					for(r = ((comp_r >= m_ratio_h_dx) ? (comp_r - m_ratio_h_dx) : 0); r <= ((comp_r <= h - m_ratio_h_dx) ? (comp_r + m_ratio_h_dx) : h); r++){
						for(int comp_index = 0; comp_index < m_grid_num[i][j][r]; comp_index++){
							int neighbor_ID = idx_table[i][j][r][comp_index];
							double r_sqnorm = (nodes[k].position - nodes[neighbor_ID].position).sqrNorm();
							if(r_sqnorm <= m_h * m_h){ nodes[k].neighbor_index.push_back(neighbor_ID); }
						}
					}
				}
			}
		}

		//step1:计算粒子的Non-Pressure加速度：
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			//nodes[k].force_vis = Vec3(0.0, 0.0, 0.0); execute in step 0.
			//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
			if(nodes[k].isWall()){ continue; }
			compute_vis_accelerate(k);
		}

		//step2:更新粒子的速度
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			if(nodes[k].isWall()){ continue; }
			nodes[k].velocity += (nodes[k].force_vis / m_mass + nodes[k].acc_ext) * m_dt;
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
				nodes[k].force_press = Vec3(0.0, 0.0, 0.0);
				nodes[k].pressure = m_KK * (pow(nodes[k].rho / rho_L, (int)m_gamma) - 1.0);
				if(nodes[k].pressure < 0.0){ nodes[k].pressure = 0.0; }
			}

			//step5:计算粒子的Pressure加速度：
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				//nodes[k].acc_press = Vec3(0.0, 0.0, 0.0); execute in step 4.
				//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
				if(!nodes[k].isFluid()){ continue; }
				compute_press_accelerate(k);
			}

			//step6:更新粒子的速度
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				if(!nodes[k].isFluid()){ continue; }
				nodes[k].velocity += nodes[k].force_press * (m_dt / m_mass);
			}

			//Check if we can break out the iteration:
			rho_error = 0.0;
#pragma omp parallel for private(k)
			for(k = 0; k < total_num; k++){
				rho_error += (nodes[k].rho - rho_L) * (nodes[k].rho - rho_L);
			}
			rho_error = sqrt(rho_error);
			if(rho_error < 0.1){ break; }
		}
		std::cerr << "rho_error: " << rho_error << std::endl;

		//step6.5:更新刚体运动
		rigidbody.update(m_dt);
		double dTheta_dt = rigidbody.swirl_velocity;
		double dTheta = dTheta_dt * m_dt;
		double cosdTheta = cos(dTheta), sindTheta = sin(dTheta);
#pragma omp parallel for private(k,i)
		for(i = 0; i < rigidbody.particle_indexes.size(); i++){
			k = rigidbody.particle_indexes[i];
			double local_x = nodes[k].position.getx() - rigidbody.position.getx();
			double local_y = nodes[k].position.gety() - rigidbody.position.gety();
			double r_norm = sqrt(local_x * local_x + local_y * local_y);
			double costheta = local_x / r_norm, sintheta = local_y / r_norm;
			//update:
			nodes[k].setPosition((costheta * cosdTheta - sintheta * sindTheta) * r_norm + rigidbody.position.getx(),
				(sintheta * cosdTheta + costheta * sindTheta) * r_norm + rigidbody.position.gety(),
				nodes[k].position.getz());
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
			if(!nodes[k].isFluid()){ continue; }
			nodes[k].position += nodes[k].velocity * m_dt;
			if(boundary_band[0] == 'P'){ //Left
				if(nodes[k].position.getx() < 0){ nodes[k].position += Vec3(x_bound, 0.0, 0.0); }
			}
			if(boundary_band[1] == 'P'){ //Right
				if(nodes[k].position.getx() > x_bound){ nodes[k].position += Vec3(-x_bound, 0.0, 0.0); }
			}
			if(boundary_band[2] == 'P'){ //Forward
				if(nodes[k].position.gety() < 0){ nodes[k].position += Vec3(0.0, y_bound, 0.0); }
			}
			if(boundary_band[3] == 'P'){ //Back
				if(nodes[k].position.gety() > y_bound){ nodes[k].position += Vec3(0.0, -y_bound, 0.0); }
			}
			if(boundary_band[4] == 'P'){ //Down
				if(nodes[k].position.getz() < 0){ nodes[k].position += Vec3(0.0, 0.0, z_bound); }
			}
			if(boundary_band[5] == 'P'){ //Up
				if(nodes[k].position.getz() > z_bound){ nodes[k].position += Vec3(0.0, 0.0, -z_bound); }
			}
		}

		Vmax = 0.0;
#pragma omp parallel for private(k)
		for(k = 0; k < total_num; k++){
			double u_maxAbsCoord = nodes[k].velocity.maxAbsCoord();
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

		fout << "H0(1,1) = " << x_bound * 100 << " A" << std::endl;
		fout << "H0(1,2) = " << 0 << " A" << std::endl;
		fout << "H0(1,3) = " << 0 << " A" << std::endl;

		fout << "H0(2,1) = " << 0 << " A" << std::endl;
		fout << "H0(2,2) = " << y_bound * 100 << " A" << std::endl;
		fout << "H0(2,3) = " << 0 << " A" << std::endl;

		fout << "H0(3,1) = " << 0 << " A" << std::endl;
		fout << "H0(3,2) = " << 0 << " A" << std::endl;
		fout << "H0(3,3) = " << z_bound * 100 << " A" << std::endl;

		fout << ".NO_VELOCITY." << std::endl;
		fout << "entry_count = " << 3 << std::endl;

		for(int k = 0; k < total_num; k++){
			fout << m_mass << std::endl; // mass
			fout << nodes[k].status << std::endl; // element type
			fout << nodes[k].position.getx() * 100 << " " << nodes[k].position.gety() * 100 << " " << nodes[k].position.getz() * 100 << std::endl; // Info of position.
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
			if(nodes[k].status == "Wall" || nodes[k].status == "RigidBody") continue; // note Wall or Rigid Body
			cnt++;
			ss << nodes[k].position.getx() << " " << nodes[k].position.gety() << " " << nodes[k].position.getz() << endl;
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
			//line3: H0(1,1) = {$x_bound * 100} A
			//line4: H0(1,2) = 0 A
			//line5: H0(1,3) = 0 A
			//line6: H0(2,1) = 0 A
			//line7: H0(2,2) = {$y_bound * 100} A
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
		std::string status;
		while(k != total_num){
			fin >> m_mass;
			fin >> status;
			nodes[k].setStatus(status);
			fin >> in_x >> in_y >> in_z;
			nodes[k].setPosition(in_x / 100.0, in_y / 100.0, in_z / 100.0);
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

		for(i = total_num - 1; i >= 0; i--){
			glColor3f(0.5f, 0.5f, 1.0f);
			if(nodes[i].isWall()){
				continue;
			}
			if(nodes[i].isRigidBody()){
				glColor3f(1.0f, 1.0f, 0.0f);
			}
			glPushMatrix();
			glTranslated(nodes[i].position.getx() / dx, nodes[i].position.gety() / dy, nodes[i].position.getz() / dz);
			glutSolidSphere(radius, 6, 6);
			glPopMatrix();
		}

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
			else if(if_dump == 3){ dump_ply_file("../output/liquid_vis0/" + ss.str() + ".ply"); }
			// in all cases dump wheel 
			// dump_wheel_ply("../output/wheel/" + ss.str() + ".ply"); // do not have to

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