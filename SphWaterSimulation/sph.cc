#include "sph.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>

#include "constants.h"

namespace SPH
{
    Sph::Sph():
        use_rigid_body(program_const::kUseRigidBody),
        if_visualize(program_const::IF_VISUALIZE),
        read_from_file(program_const::READ_FROM_FILES),
        dump_file_interval(program_const::kDumpFileInterval),
        dh_ratio(sph_const::kdhRatio),
        lambda(sph_const::kLambda),
        stiffness(sph_const::kStiffness),
        mass_liquid(sph_const::kMassLiquid),
        gamma(sph_const::kGamma),
        viscosity(sph_const::kViscosity){

        cur_step = 0;
        nodes_num = 0;
        init_pressure = 0.;
        max_velocity = 0.;
        dx = x_bound / l, dy = y_bound / m, dz = z_bound / h;
        dh = dh_ratio * dx;
        rho_liquid = mass_liquid / dx / dy / dz;
        dt = lambda * dh / sqrt(stiffness);
        KK = stiffness * rho_liquid / gamma;
        memset(grid_num, 0, sizeof(grid_num));
        rigidbody = Wheel(sph_const::kRhoRigidbody, sph_const::kRigidBodyRadiusOutSize,
            sph_const::kRigidBodyRadiusInSize, sph_const::kRigidBodyHeight,
            sph_const::kRigidBodyCenter, sph_const::kRigidBodyLeafNum);
        initNodes();
        setKernels();
    }

    void Sph::initNodes(){
        // Set the bounding box as "Wall" nom_hdes
        setWallNodes();
        // Set the wheel as "RigidBody" nodes, as well as the meshes for rendering.
        if(use_rigid_body){
            setRigidBodyMesh();
            setRigidBodyNodes();
        }
        // Set the water as "Liquid" nodes
        setLiquidNodes();
        std::cout << "The total number of particles: " << nodes_num << std::endl;
        prev_velocity = std::vector<Vec3>(nodes_num, Vec3(0.0, 0.0, 0.0));
    }

    void Sph::setKernels(){
        double m_Kpoly6;
        double a_Kpoly6;
        double vis_Kpoly6;
        m_Kpoly6 = 315 / ((64 * math_const::PI) * pow(dh, 9));
        a_Kpoly6 = 45 / (math_const::PI * pow(dh, 6));
        vis_Kpoly6 = 15.0 / (2.0 * math_const::PI * pow(dh, 3));
        kernels = Kernels(dh, m_Kpoly6, a_Kpoly6, vis_Kpoly6);
    }

    void Sph::step(){
        if(read_from_file){
            readFileAndShow();
        }
        else{
            if(if_visualize)
                draw();
            compute();
            std::cout << "step: " << cur_step << "\terror: " << getError();
            if(cur_step++ % dump_file_interval == 0){
                dumpFiles();
            }
        }
    }

    void Sph::readFileAndShow(){
        if(program_const::IF_VISUALIZE){
            draw();
        }
        if(cur_step++ % dump_file_interval == 0){
            std::string file_path = getFilePath("readLiquidFromCfg");
            if(!readLiquidFromCfg(file_path)){
                cur_step = 0;
            }
        }
    }

    void Sph::compute(){
        int i = 0, j = 0, r = 0, k = 0;

        // Step 0: update tables
#pragma omp parallel for private(i,j,r)
        for(i = 0; i <= l; i++){
            for(j = 0; j <= m; j++){
                for(r = 0; r <= h; r++){
                    idx_table[i][j][r].clear();
                    grid_num[i][j][r] = 0;
                }
            }
        }
        for(k = 0; k < nodes_num; k++){
            i = (int)(nodes[k].position.getx() / dx);
            j = (int)(nodes[k].position.gety() / dy);
            r = (int)(nodes[k].position.getz() / dz);
            if(i < 0 || i > l || j < 0 || j > m || r < 0 || r > h){
                continue;
            }
            idx_table[i][j][r].push_back(k);
            grid_num[i][j][r] += 1;
        }
        for(k = 0; k < nodes_num; k++){
            nodes[k].force_vis = Vec3(0.0, 0.0, 0.0);
            nodes[k].neighbor_index.clear();
            int comp_i = (int)(nodes[k].position.getx() / dx);
            int comp_j = (int)(nodes[k].position.gety() / dy);
            int comp_r = (int)(nodes[k].position.getz() / dz);
            for(i = ((comp_i >= dh_ratio) ? (comp_i - dh_ratio) : 0); i <= ((comp_i <= l - dh_ratio) ? (comp_i + dh_ratio) : l); i++){
                for(j = ((comp_j >= dh_ratio) ? (comp_j - dh_ratio) : 0); j <= ((comp_j <= m - dh_ratio) ? (comp_j + dh_ratio) : m); j++){
                    for(r = ((comp_r >= dh_ratio) ? (comp_r - dh_ratio) : 0); r <= ((comp_r <= h - dh_ratio) ? (comp_r + dh_ratio) : h); r++){
                        for(int comp_index = 0; comp_index < grid_num[i][j][r]; comp_index++){
                            int neighbor_ID = idx_table[i][j][r][comp_index];
                            double r_sqnorm = (nodes[k].position - nodes[neighbor_ID].position).sqrNorm();
                            if(r_sqnorm <= dh * dh){ nodes[k].neighbor_index.push_back(neighbor_ID); }
                        }
                    }
                }
            }
        }

        // Step 1: Compute non-pressure accelerate for particles
#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            //nodes[k].force_vis = Vec3(0.0, 0.0, 0.0); execute in step 0.
            if(nodes[k].isWall()){ continue; }
            computeVisAccerlerate(k);
        }

        // Step 2: Update particle velocity
#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            if(nodes[k].isWall()){ continue; }
            nodes[k].velocity += (nodes[k].force_vis / mass_liquid + nodes[k].acc_ext) * dt;
        }

        // Step 3: Compute density with kernel functions
        int iter_num = 0;
        double rho_error = 0.0;
        while(++iter_num <= 10){
#pragma omp parallel for private(k)
            for(k = 0; k < nodes_num; k++){
                computeDensity(k);
            }

            // Step 4: Compute pressures for particles
#pragma omp parallel for private(k)
            for(k = 0; k < nodes_num; k++){
                nodes[k].force_press = Vec3(0.0, 0.0, 0.0);
                nodes[k].pressure = KK * (pow(nodes[k].rho / rho_liquid, (int)gamma) - 1.0);
                if(nodes[k].pressure < 0.0){ nodes[k].pressure = 0.0; }
            }
            // Step 5: Compute pressure accelerate for particles
#pragma omp parallel for private(k)
            for(k = 0; k < nodes_num; k++){
                // Compute only pressure accelerate for only liquid,
                // since that of rigidbody can got from the reaction force,
                // and that of boundary won't move
                if(!nodes[k].isLiquid()){ continue; }
                computePressureAccelerate(k);
            }

            // Step 6: Update particle velocity
#pragma omp parallel for private(k)
            for(k = 0; k < nodes_num; k++){
                if(!nodes[k].isLiquid()){ continue; }
                nodes[k].velocity += nodes[k].force_press * (dt / mass_liquid);
            }

            // Check if we can break the iteration
            rho_error = 0.0;
#pragma omp parallel for private(k)
            for(k = 0; k < nodes_num; k++){
                rho_error += (nodes[k].rho - rho_liquid) * (nodes[k].rho - rho_liquid);
            }
            rho_error = sqrt(rho_error);
            if(rho_error < 0.1){ break; }
        }
        std::cout << "\trho_error: " << rho_error << std::endl;

        // Step 7: Update rigidbody
        if(use_rigid_body){
            updateRigidBody();
        }

        // Step 8: Update particle positions
#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            if(!nodes[k].isLiquid()){ continue; }
            nodes[k].position += nodes[k].velocity * dt;
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

        max_velocity = 0.0;
#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            double u_maxAbsCoord = nodes[k].velocity.maxAbsCoord();
            max_velocity = max_velocity > u_maxAbsCoord ? max_velocity : u_maxAbsCoord;
        }
        // Update dt
        dt = lambda * dh / std::max(max_velocity, sqrt(stiffness));
        dt = lambda * dh / sqrt(stiffness);

        return;
    }

    void Sph::computeVisAccerlerate(int k){
        double Vb = 0.0;
        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];
            // Fix density for moving rigidbody and boundary, Akinci et al. (2012)
            if(!nodes[neighbor_ID].isLiquid()){
                Vb += kernels.Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
            }
        }
        Vb = 1.0 / Vb;

        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];
            if(neighbor_ID == k){ continue; }

            double r_sqnorm = (nodes[k].position - nodes[neighbor_ID].position).sqrNorm();
            if(r_sqnorm > dh * dh){ continue; }

            Vec3 Gradient_Wij = kernels.Gradient_Muller03Kernel_Pressure(nodes[k].position, nodes[neighbor_ID].position);
            double u_dot_r = (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(nodes[k].position - nodes[neighbor_ID].position);
            if(nodes[neighbor_ID].isLiquid()){
                //fluid - fluid
                nodes[k].force_vis += Gradient_Wij * ((viscosity) * 2.0 * (3 + 2.0) * (mass_liquid * mass_liquid / nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
            }
            else{
                //fluid - obstacle
                double mass_rigid = mass_liquid;
                Vec3 temp_force = Gradient_Wij * ((viscosity) * 2.0 * (3 + 2.0) * (mass_rigid * mass_liquid / nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
                nodes[k].force_vis += temp_force;
                nodes[neighbor_ID].force_vis += -temp_force;
            }
        }
    }

    void Sph::computePressureAccelerate(int k){
        if(nodes[k].pressure < 0.0){ return; }
        double Vb = 0.0;
        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];
            // Fix density for moving rigidbody and boundary, Akinci et al. (2012)
            if(!nodes[neighbor_ID].isLiquid()){
                Vb += kernels.Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
            }
        }
        Vb = 1.0 / Vb;

        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];

            if(neighbor_ID == k){ continue; }
            if((nodes[k].position - nodes[neighbor_ID].position).sqrNorm() > dh * dh){ continue; }

            Vec3 Gradient_Wij = kernels.Gradient_Muller03Kernel_Pressure(nodes[k].position, nodes[neighbor_ID].position);
            if(nodes[neighbor_ID].isWall()){
                //fluid - Obstacle
                double mass_rigid = mass_liquid;
                Vec3 temp_force = Gradient_Wij * (-2.0 * mass_liquid * mass_rigid * nodes[k].pressure / nodes[k].rho / nodes[k].rho);
                nodes[k].force_press += temp_force;
                nodes[neighbor_ID].force_press += -temp_force;
            }
            else{
                //fluid - fluid
                nodes[k].force_press += Gradient_Wij * (-mass_liquid * mass_liquid * \
                    (nodes[k].pressure / nodes[k].rho / nodes[k].rho + \
                        nodes[neighbor_ID].pressure / nodes[neighbor_ID].rho / nodes[neighbor_ID].rho));
            }
        }
    }

    void Sph::computeDensity(int k){
        if(nodes[k].isWall()){ return; }
        nodes[k].rho = 0.0;
        double Vb = 0.0;
        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];
            // Fix density for moving rigidbody and boundary, Akinci et al. (2012)
            if(!nodes[neighbor_ID].isLiquid()){
                Vb += kernels.Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
            }
        }
        Vb = 1.0 / Vb;
        for(int i = 0; i < nodes[k].neighbor_index.size(); i++){
            int neighbor_ID = nodes[k].neighbor_index[i];
            if((nodes[k].position - nodes[neighbor_ID].position).sqrNorm() > dh * dh){ continue; }

            Vec3 Grad_W = kernels.Gradient_Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
            if(!nodes[neighbor_ID].isLiquid()){
                double mass_rigid = mass_liquid;
                nodes[k].rho += mass_rigid * kernels.Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
                nodes[k].rho += mass_rigid * dt * (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(Grad_W);
            }
            else{
                nodes[k].rho += mass_liquid * kernels.Muller03Kernel_Basic(nodes[k].position, nodes[neighbor_ID].position);
                nodes[k].rho += mass_liquid * dt * (nodes[k].velocity - nodes[neighbor_ID].velocity).dot(Grad_W);
            }
        }
    }

    void Sph::updateRigidBody(){
        rigidbody.update(dt);
        double dTheta_dt = rigidbody.swirl_velocity;
        double dTheta = dTheta_dt * dt;
        double cosdTheta = cos(dTheta), sindTheta = sin(dTheta);
        int k, i;
#pragma omp parallel for private(k,i)
        for(i = 0; i < rigidbody.particle_indexes.size(); i++){
            k = rigidbody.particle_indexes[i];
            double local_x = nodes[k].position.getx() - rigidbody.position.getx();
            double local_y = nodes[k].position.gety() - rigidbody.position.gety();
            double r_norm = sqrt(local_x * local_x + local_y * local_y);
            double costheta = local_x / r_norm, sintheta = local_y / r_norm;
            nodes[k].setPosition((costheta * cosdTheta - sintheta * sindTheta) * r_norm + rigidbody.position.getx(),
                (sintheta * cosdTheta + costheta * sindTheta) * r_norm + rigidbody.position.gety(),
                nodes[k].position.getz());
        }
        Mat3 RotationMat(std::vector<Vec3>{Vec3(cosdTheta, -sindTheta, 0.0),
            Vec3(sindTheta, cosdTheta, 0.0),
            Vec3(0.0, 0.0, 1.0)});
        for(i = 0; i < localcoord_rigidbody.size(); i++){
            localcoord_rigidbody[i] = RotationMat.multiply(localcoord_rigidbody[i]);
        }
    }

    // Register a node of different kinds to 'nodes' as well as 'idx_table' during initNodes().
    // If the node is rigidbody
    void Sph::addNode(int i, int j, int r, std::string status){
        // Rigidbody and wall nodes takes the first priority.
        if(grid_num[i][j][r] == 1){
            if(status != "Liquid"){
                if(!nodes[idx_table[i][j][r][0]].isRigidBody() && status == "RigidBody"){
                    rigidbody.particle_indexes.push_back(idx_table[i][j][r][0]);
                }
                nodes[idx_table[i][j][r][0]].setStatus(status);
            }
            return;
        }
        Particle temp_node(rho_liquid, init_pressure, g_gravity);
        // Register the index of particle to 'idx_table'.
        idx_table[i][j][r].push_back(nodes_num++);
        temp_node.setPosition(i * dx, j * dy, r * dz);
        temp_node.setStatus(status);
        nodes.push_back(temp_node);//Push new node into node list.
        if(temp_node.isRigidBody()){ rigidbody.particle_indexes.push_back(nodes.size() - 1); }
        grid_num[i][j][r] = 1;
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
        translate_rigidbody = rigidbody.position; translate_rigidbody.setz(0.0);
        int i = 0, j = 0, r = 0;
        for(i = 0; i <= l; i++){
            for(j = 0; j <= m; j++){
                for(r = 0; r <= h; r++){
                    if(r * dz > 2 * dz + rigidbody.m_wheel_height){ continue; }  // 2*dy due to the boundary layer.
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

    void Sph::setLiquidNodes(){
        for(int i = 0; i <= l; i++){
            for(int j = 0; j < m / 3; j++){
                for(int r = 0; r < 3 * h / 4; r++){
                    addNode(i, j, r, "Liquid");
                }
            }
        }
    }

    void Sph::setRigidBodyMesh(){
        double down_height = 2 * dz, up_height = rigidbody.m_wheel_height + 2 * dz;
        double r_in = rigidbody.m_wheel_radius_insize, r_out = rigidbody.m_wheel_radius_outsize;
        localcoord_rigidbody.push_back(Vec3(0.0f, 0.0f, down_height));  // 2*dy is due to the boundary layer
        localcoord_rigidbody.push_back(Vec3(0.0f, 0.0f, up_height));  // 2*dy is due to the boundary layer
        double dtheta = math_const::PI / rigidbody.m_leafnum;
        for(int k = 1; k <= rigidbody.m_leafnum; k++){
            // Step 1: Add 8 new vertex into the list
            int idx_8[8] = { 8 * k - 6, 8 * k - 5, 8 * k - 4, 8 * k - 3, 8 * k - 2, 8 * k - 1, 8 * k, 8 * k + 1 };
            double cos_thetaminor = cos((2 * k - 2) * dtheta), sin_thetaminor = sin((2 * k - 2) * dtheta);
            double cos_thetabigger = cos((2 * k - 1) * dtheta), sin_thetabigger = sin((2 * k - 1) * dtheta);

            localcoord_rigidbody.push_back(Vec3(r_in * cos_thetaminor, r_in * sin_thetaminor, down_height));  // idx_8[0]
            localcoord_rigidbody.push_back(Vec3(r_out * cos_thetaminor, r_out * sin_thetaminor, down_height));  // idx_8[1]
            localcoord_rigidbody.push_back(Vec3(r_out * cos_thetabigger, r_out * sin_thetabigger, down_height));  // idx_8[2]
            localcoord_rigidbody.push_back(Vec3(r_in * cos_thetabigger, r_in * sin_thetabigger, down_height));  // idx_8[3]

            localcoord_rigidbody.push_back(Vec3(r_in * cos_thetaminor, r_in * sin_thetaminor, up_height));  // idx_8[4]
            localcoord_rigidbody.push_back(Vec3(r_out * cos_thetaminor, r_out * sin_thetaminor, up_height));  // idx_8[5]
            localcoord_rigidbody.push_back(Vec3(r_out * cos_thetabigger, r_out * sin_thetabigger, up_height));  // idx_8[6]
            localcoord_rigidbody.push_back(Vec3(r_in * cos_thetabigger, r_in * sin_thetabigger, up_height));  // idx_8[7]

            // Step 2: Create faces_rigidbody: 4 triangles + 4 Quads = 12 * triangles
            faces_rigidbody.push_back(std::vector<int>{0, idx_8[2], idx_8[1]});  // Trigngle #1
            faces_rigidbody.push_back(std::vector<int>{1, idx_8[5], idx_8[6]});  // Trigngle #2

            faces_rigidbody.push_back(std::vector<int>{idx_8[1], idx_8[2], idx_8[6]});  // Quad #0.5
            faces_rigidbody.push_back(std::vector<int>{idx_8[1], idx_8[6], idx_8[5]});  // Quad #1
            faces_rigidbody.push_back(std::vector<int>{idx_8[2], idx_8[3], idx_8[7]});  // Quad #1.5
            faces_rigidbody.push_back(std::vector<int>{idx_8[2], idx_8[7], idx_8[6]});  // Quad #2
            faces_rigidbody.push_back(std::vector<int>{idx_8[0], idx_8[1], idx_8[5]});  // Quad #2.5
            faces_rigidbody.push_back(std::vector<int>{idx_8[0], idx_8[5], idx_8[4]});  // Quad #3

            int idx_next0 = idx_8[0] + 8, idx_next4 = idx_8[4] + 8;
            if(idx_next4 >= 2 + 8 * rigidbody.m_leafnum){
                idx_next0 -= rigidbody.m_leafnum * 8;
                idx_next4 -= rigidbody.m_leafnum * 8;
            }
            faces_rigidbody.push_back(std::vector<int>{0, idx_next0, idx_8[3]});  // Trigngle #3
            faces_rigidbody.push_back(std::vector<int>{1, idx_8[7], idx_next4});  // Trigngle #4

            faces_rigidbody.push_back(std::vector<int>{idx_8[3], idx_next0, idx_next4}); //Quad #3.5
            faces_rigidbody.push_back(std::vector<int>{idx_8[3], idx_next4, idx_8[7]}); //Quad #4
        }
    }

    void Sph::dumpFiles(){
        std::string file_path;
        if(program_const::kDumpLiquidAsCfg){
            file_path = getFilePath("dumpLiquidAsCfg");
            std::cout << "dump " << file_path << std::endl;
            dumpLiquidAsCfg(file_path);
        }
        if(program_const::kDumpLiquidAsPly){
            file_path = getFilePath("dumpLiquidAsPly");
            std::cout << "dump " << file_path << std::endl;
            dumpLiquidAsPly(file_path);
        }
        if(program_const::kDumpRigidBody){
            file_path = getFilePath("dumpRigidBoty");
            std::cout << "dump " << file_path << std::endl;
            dumpRigidBody(file_path);
        }
    }
    std::string Sph::getFilePath(std::string command){
        std::stringstream ss;
        if(command == "dumpLiquidAsPly"){
            ss << program_const::PLY_SUBFOLDER;
            ss << std::setfill('0') << std::setw(4) << cur_step / dump_file_interval;
            ss << ".ply";
        }
        else if(command == "dumpRigidBody"){
            ss << program_const::RIGIDBODY_SUBFOLDER;
            ss << std::setfill('0') << std::setw(4) << cur_step / dump_file_interval;
            ss << ".ply";
        }
        else if(command == "dumpLiquidAsCfg" || command == "readLiquidFromCfg"){
            ss << program_const::CFG_SUBFOLDER;
            ss << std::setfill('0') << std::setw(4) << cur_step / dump_file_interval;
            ss << ".cfg";
        }
        else{
            std::cerr << "[getFilePath] called with known command: " << command << std::endl;
        }
        return ss.str();
    }

    // This function can be used to dump a file which can be visualized by software OVITO.
    void Sph::dumpLiquidAsCfg(std::string file_path){
        std::ofstream fout(file_path);
        if(!fout){
            std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
            exit(-1);
        }

        fout << "Number of particles = " << nodes_num << std::endl;
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

        for(int k = 0; k < nodes_num; k++){
            fout << mass_liquid << std::endl; // mass
            fout << nodes[k].status << std::endl; // element type
            fout << nodes[k].position.getx() * 100 << " " << nodes[k].position.gety() * 100 << " " << nodes[k].position.getz() * 100 << std::endl; // Info of position.
        }

        fout.close();

        return;
    }

    void Sph::dumpRigidBody(std::string file_path){
        if(!use_rigid_body){
            return;
        }
        std::ofstream fout(file_path);
        if(!fout){
            std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
            exit(-1);
        }
        fout << "ply" << std::endl;
        fout << "format ascii 1.0" << std::endl;
        fout << "comment made by anonymous" << std::endl;
        fout << "comment this file is a wheel" << std::endl;
        fout << "element vertex " << localcoord_rigidbody.size() << std::endl;
        fout << "property float32 x" << std::endl;
        fout << "property float32 y" << std::endl;
        fout << "property float32 z" << std::endl;
        fout << "element face " << faces_rigidbody.size() << std::endl;
        fout << "property list uint8 int32 vertex_index" << std::endl;
        fout << "end_header" << std::endl;

        // Dump vertex coordinate into file
        for(int k = 0; k < localcoord_rigidbody.size(); k++){
            Vec3 single_vertex_worldcoord = localcoord_rigidbody[k] + translate_rigidbody;
            fout << single_vertex_worldcoord.getx() << " " << single_vertex_worldcoord.gety() << " " << \
                single_vertex_worldcoord.getz() << std::endl;
        }
        // Dump the vertex index of faces_rigidbody into file
        for(int k = 0; k < faces_rigidbody.size(); k++){
            fout << "3 ";
            for(int rr = 0; rr < faces_rigidbody[k].size(); rr++){
                fout << faces_rigidbody[k][rr] << " ";
            }
            fout << std::endl;
        }
    }

    void Sph::dumpLiquidAsPly(std::string file_path){
        std::ofstream fout(file_path);
        if(!fout){
            std::cerr << "[dumpLiquidAsPly(" << file_path << ")] fails!" << std::endl;
            exit(-1);
        }
        std::stringstream ss;
        int vertexCnt = 0;
        for(int i = 0; i < nodes_num; i++){
            if(nodes[i].isLiquid()){
                vertexCnt++;
                const Vec3& position = nodes[i].position;
                ss << position.getx() << " " << position.gety() << " " << position.getz() << std::endl;
            }
        }
        fout << "ply" << std::endl;
        fout << "format ascii 1.0" << std::endl;
        fout << "element vertex " << vertexCnt << std::endl;
        fout << "property float x" << std::endl;
        fout << "property float y" << std::endl;
        fout << "property float z" << std::endl;
        fout << "end_header" << std::endl;
        fout << ss.str();
    }

    // Returns true if reading success
    bool Sph::readLiquidFromCfg(std::string file_path){
        std::ifstream fin(file_path);
        std::string get_str;
        if(!fin){
            std::cerr << "[readLiquidFromCfg(" << file_path << ")] fails!" << std::endl;
            return false;
        }
        for(int k = 0; k < 13; k++){
            getline(fin, get_str);
        }
        int k = 0;
        double in_x, in_y, in_z;
        std::string status;
        while(k != nodes_num){
            fin >> mass_liquid;
            fin >> status;
            nodes[k].setStatus(status);
            fin >> in_x >> in_y >> in_z;
            nodes[k].setPosition(in_x / 100.0, in_y / 100.0, in_z / 100.0);
            k++;
        }
        std::cout << "cur_step: " << cur_step << ", nodes_num: " << nodes_num << std::endl;
        return true;
    }

    void Sph::draw(){
        assert(if_visualize && "Calling draw() with IF_VISUALIZE == false!");
        drawBoundingBox();
        drawLiquid();
        drawRigidBody();
        return;
    }

    void Sph::drawBoundingBox(){
        glColor3f(0.0, 1.0, 0.0);
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
    }

    void Sph::drawRigidBody(){
        glColor3f(1.0f, 1.0f, 0.0f);
        double size = visual_const::kParticleRenderingSize;
        for(int i = nodes_num - 1; i >= 0; i--){
            if(nodes[i].isRigidBody()){
                glPushMatrix();
                glTranslated(nodes[i].position.getx() / dx,
                    nodes[i].position.gety() / dy, nodes[i].position.getz() / dz);
                glutSolidSphere(size, 6, 6);
                glPopMatrix();
            }
        }
    }

    void Sph::drawLiquid(){
        glColor3f(0.5f, 0.5f, 1.0f);
        double size = visual_const::kParticleRenderingSize;
        for(int i = nodes_num - 1; i >= 0; i--){
            if(nodes[i].isLiquid()){
                glPushMatrix();
                glTranslated(nodes[i].position.getx() / dx,
                    nodes[i].position.gety() / dy, nodes[i].position.getz() / dz);
                glutSolidSphere(size, 6, 6);
                glPopMatrix();
            }
        }
    }

    void Sph::recordVelocity(){
        int k = 0;
#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            prev_velocity[k] = nodes[k].velocity;
        }
    }

    double Sph::getError(){
        int k = 0;
        double final_error = 0.0;

#pragma omp parallel for private(k)
        for(k = 0; k < nodes_num; k++){
            final_error += (prev_velocity[k] - nodes[k].velocity).sqrNorm();
        }
        return sqrt(final_error);
    }

    int Sph::getStep(){
        return cur_step;
    }
}