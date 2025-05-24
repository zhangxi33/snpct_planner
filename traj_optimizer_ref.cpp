#include "plan_manage/traj_optimizer.h"
// using namespace std;

namespace plan_manage
{
    
//
 double  PolyTrajOptimizer::costFunctionCallback_ref_RLC(void *func_data, const Eigen::VectorXd &x_l, Eigen::VectorXd &grad_x_l){
      double total_smcost = 0.0, total_timecost = 0.0, penalty_cost = 0.0, totoal_refcost = 0.0;
       PolyTrajOptimizer *opt = reinterpret_cast<PolyTrajOptimizer *>(func_data);

        int piece_num_ = opt->piece_num_container[0];

        grad_x_l.resize(piece_num_-1+opt->trajnum);
        grad_x_l.setZero();

        int variable_num_ = (piece_num_-1)*2;
        variable_num_ += opt->trajnum;
        Eigen::VectorXd x;
        Eigen::VectorXd grad;


        x.resize(variable_num_);
        grad.resize(variable_num_);

    

     for(int i = 0; i < x_l.size()-1; i++){
         x[i*2] = opt->initInnerPts_ref_container[0](1,i)-x_l(i)*sin(opt->initInnerPts_ref_container[0](3,i));
         x[i*2+1] = opt->initInnerPts_ref_container[0](2,i)+x_l(i)*cos(opt->initInnerPts_ref_container[0](3,i));


        //  std::cout<<"x_l(i): "<<x_l(i)
        //          <<" ; "<<x[i*2]
        //          <<" ; "<<x[i*2+1]
        //            <<std::endl;
     }
     x[x.size()-1]=x_l[x_l.size()-1];
        
        

       std::vector<Eigen::Map<const Eigen::MatrixXd>> P_container;
       std::vector<Eigen::Map<Eigen::MatrixXd>> gradP_container;;
       int offset = 0;

// //
       for(int trajid = 0; trajid < opt->trajnum; trajid++){
           Eigen::Map<const Eigen::MatrixXd> P(x.data()+offset, 2, opt->piece_num_container[trajid] - 1);
           Eigen::Map<Eigen::MatrixXd> gradP(grad.data()+offset, 2, opt->piece_num_container[trajid] - 1);
           offset += 2 * (opt->piece_num_container[trajid] - 1);
           gradP.setZero();
           P_container.push_back(P);
           gradP_container.push_back(gradP);
       }
       Eigen::Map<const Eigen::VectorXd> t(x.data()+offset, opt->trajnum);
       Eigen::Map<Eigen::VectorXd> gradt(grad.data()+offset, opt->trajnum);
       offset += opt -> trajnum;
       Eigen::VectorXd T(opt->trajnum);
       Eigen::VectorXd gradT(opt->trajnum); 
       gradT.setZero();
       opt->VirtualT2RealT(t, T);
//        std::vector<double> trajtimes; 
//        trajtimes.push_back(0.0);
// //        //T(i) is sum time of i-segment traj
//        for(int trajid = 0; trajid < opt->trajnum; trajid++){
//            trajtimes.push_back(T[trajid]);
//        }
//        // Ini/Fin Gear Pos
//        std::vector<Eigen::Map<const Eigen::MatrixXd>> Gear_container;
//        std::vector<Eigen::Map<Eigen::MatrixXd>> gradGear_container;
//        for(int trajid = 0; trajid < opt->trajnum - 1; trajid++){
//            Eigen::Map<const Eigen::MatrixXd> Gear(x.data()+offset, 2, 1);
//            Eigen::Map<Eigen::MatrixXd>gradGear(grad.data()+offset, 2, 1);
//            offset += 2;
//            gradGear.setZero();
//            Gear_container.push_back(Gear);
//            gradGear_container.push_back(gradGear);
//        }
// //        //
//        Eigen::Map<const Eigen::VectorXd> Angles(x.data()+offset, opt->trajnum-1);
//        Eigen::Map<Eigen::VectorXd>gradAngles(grad.data()+offset, opt->trajnum-1);
//        gradAngles.setZero();
//
//        /*for relax end v?*/
//        // offset += opt->trajnum-1;
//        // Eigen::Map<const Eigen::MatrixXd> endV (x.data()+offset, 2, 1);
//        // Eigen::Map<Eigen::MatrixXd> gradendV (grad.data()+offset, 2, 1);
//        // gradendV.setZero();
//
//
//
//
//
//
//
       // Eigen::VectorXd gradt;
       for(int trajid = 0; trajid < opt->trajnum; trajid++){
           double smoo_cost;
           Eigen::VectorXd obs_surround_feas_qvar_costs(4);
           obs_surround_feas_qvar_costs.setZero();

           Eigen::MatrixXd IniS,FinS;
           IniS = opt->iniState_container[trajid];
           FinS = opt->finState_container[trajid];

        //    if(trajid > 0){
        //        double theta = Angles[trajid-1];
        //        IniS.col(0) = Gear_container[trajid-1];
            //    IniS.col(1) = Eigen::Vector2d(-opt->non_sinv*cos(theta), -opt->non_sinv*sin(theta));
        //    }
        //    if(trajid < opt->trajnum-1){
        //        double theta = Angles[trajid];
        //        FinS.col(0) = Gear_container[trajid];
        //        FinS.col(1) = Eigen::Vector2d(opt->non_sinv*cos(theta), opt->non_sinv*sin(theta));
        //    }
// //            //relax end v?
// //            // if(trajid == opt->trajnum - 1){
// //            //   FinS.col(1) = endV;
// //            // }
// //
           opt->jerkOpt_container[trajid].generate(P_container[trajid],T[trajid] / opt->piece_num_container[trajid],IniS,FinS);
           opt->jerkOpt_container[trajid].initSmGradCost(); // Smoothness cost
           smoo_cost = opt->jerkOpt_container[trajid].getTrajJerkCost();
           opt->addPVAGradCost2CT( obs_surround_feas_qvar_costs, trajid, 0.0); // Time int cost
           //Get gradT gradC
           total_smcost += smoo_cost;
           penalty_cost +=  obs_surround_feas_qvar_costs.sum();
// //            // std::cout<<"Trajid: "<<trajid<<" penalty: "<<obs_surround_feas_qvar_costs.transpose()<<std::endl;
       }
// //
// //
       for(int trajid = 0; trajid < opt->trajnum; trajid++){
           double time_cost = 0.0;
        //    Eigen::Matrix<double,2,3> gradIni, gradFin;
           opt->jerkOpt_container[trajid].calGrads_PT(); // gdt gdp gdhead gdtail
           //waypoint
           gradP_container[trajid] = opt->jerkOpt_container[trajid].get_gdP();
           //init Fin
        //    gradIni = opt->jerkOpt_container[trajid].get_gdHead();
        //    gradFin = opt->jerkOpt_container[trajid].get_gdTail();
        //    if(opt->GearOpt){
        //        if(trajid > 0){
        //            double theta = Angles[trajid-1];
        //            gradGear_container[trajid-1] += gradIni.col(0);
        //            gradAngles[trajid-1] += gradIni.col(1).transpose() * Eigen::Vector2d(opt->non_sinv * sin(theta), -opt->non_sinv*cos(theta));
        //        }
        //        if(trajid < opt->trajnum-1){
        //            double theta = Angles[trajid];
        //            gradGear_container[trajid] += gradFin.col(0);
        //            gradAngles[trajid] += gradFin.col(1).transpose() * Eigen::Vector2d(-opt->non_sinv * sin(theta), opt->non_sinv*cos(theta));
        //        }
        //    }

           // if(trajid ==  opt->trajnum-1){
           //   gradendV += gradFin.col(1);
           // }

           opt->VirtualTGradCost(T[trajid],t[trajid],opt->jerkOpt_container[trajid].get_gdT() / opt->piece_num_container[trajid],gradt[trajid],time_cost);

           // gradt[trajid] = 0.0;

           total_timecost += time_cost;
           // std::cout<<"gradp: \n"<<gradP_container[trajid] <<std::endl;
           // std::cout<<"gradeigen p: \n"<<opt->jerkOpt_container[trajid].get_gdP()<<std::endl;
       }
       
    //    std::cout<<"opt->wei_ref_: " <<opt->wei_ref_<<std::endl;
        for(int i = 0; i<opt->piece_num_container[0]-1; i++){
            totoal_refcost += x_l[i]*x_l[i];
            grad_x_l[i]=-sin(opt->initInnerPts_ref_container[0](3,i))*gradP_container[0](0,i)+
            cos(opt->initInnerPts_ref_container[0](3,i))*gradP_container[0](1,i);
            grad_x_l[i] += opt->wei_ref_*2*x_l[i];
        }
        grad_x_l[opt->piece_num_container[0]-2+1]=gradt[0];

        totoal_refcost = totoal_refcost *opt->wei_ref_;
// std::cout<<"totoal_refcost: " << totoal_refcost <<std::endl;
       opt->iter_num_ += 1;
       // std::cout << "grad angle : "<< gradAngles.transpose() << "\n";
       // std::cout << "angle : "<< Angles.transpose() << "\n";
       // std::cout <<"endV: " <<endV.transpose()<<"\n";

       return total_smcost + total_timecost + penalty_cost 
       + totoal_refcost
       ;
 }

 double  PolyTrajOptimizer::costFunctionCallback_ref_RLC_UN(void *func_data, const Eigen::VectorXd &x_l, Eigen::VectorXd &grad_x_l){
    // std::cout<<" enter costFunctionCallback_ref_RLC_UN"<<std::endl;
    double total_smcost = 0.0, total_timecost = 0.0, penalty_cost = 0.0, totoal_refcost = 0.0;
     PolyTrajOptimizer *opt = reinterpret_cast<PolyTrajOptimizer *>(func_data);

      int piece_num_ = opt->piece_num_container[0];

      grad_x_l.resize(piece_num_-1 + piece_num_);
      grad_x_l.setZero();

      int variable_num_ = (piece_num_-1)*2;
      variable_num_ += piece_num_;
      Eigen::VectorXd x;
      Eigen::VectorXd grad;


      x.resize(variable_num_);
      grad.resize(variable_num_);

  

        for(int i = 0; i < piece_num_-1; i++){
            x[i*2] = opt->initInnerPts_ref_container[0](1,i)-x_l(i)*sin(opt->initInnerPts_ref_container[0](3,i));
            x[i*2+1] = opt->initInnerPts_ref_container[0](2,i)+x_l(i)*cos(opt->initInnerPts_ref_container[0](3,i));


            //  std::cout<<"x_l(i): "<<x_l(i)
            //          <<" ; "<<x[i*2]
            //          <<" ; "<<x[i*2+1]
            //            <<std::endl;
        }
        // std::cout<<" 2"<<std::endl;

        x.tail(piece_num_) = x_l.tail(piece_num_);      
      

     std::vector<Eigen::Map<const Eigen::MatrixXd>> P_container;
     std::vector<Eigen::Map<Eigen::MatrixXd>> gradP_container;;
     int offset = 0;
    //  std::cout<<" 3"<<std::endl;

// //
     for(int trajid = 0; trajid < opt->trajnum; trajid++){
         Eigen::Map<const Eigen::MatrixXd> P(x.data()+offset, 2, opt->piece_num_container[trajid] - 1);
         Eigen::Map<Eigen::MatrixXd> gradP(grad.data()+offset, 2, opt->piece_num_container[trajid] - 1);
         offset += 2 * (opt->piece_num_container[trajid] - 1);
         gradP.setZero();
         P_container.push_back(P);
         gradP_container.push_back(gradP);
     }
     Eigen::Map<const Eigen::VectorXd> t(x.data()+offset, piece_num_);
     Eigen::Map<Eigen::VectorXd> gradt(grad.data()+offset, piece_num_);
     offset += piece_num_;
     Eigen::VectorXd T(piece_num_);

    //  Eigen::VectorXd gradT(piece_num_); 
    //  gradT.setZero();
     opt->VirtualT2RealT(t, T);
    //  std::vector<double> trajtimes; 
    //  trajtimes.push_back(0.0);
//        //T(i) is sum time of i-segment traj
    //  for(int trajid = 0; trajid < opt->trajnum; trajid++){
        //  trajtimes.push_back(T[trajid]);
    //  }
//        // Ini/Fin Gear Pos
//      std::vector<Eigen::Map<const Eigen::MatrixXd>> Gear_container;
//      std::vector<Eigen::Map<Eigen::MatrixXd>> gradGear_container;
//      for(int trajid = 0; trajid < opt->trajnum - 1; trajid++){
//          Eigen::Map<const Eigen::MatrixXd> Gear(x.data()+offset, 2, 1);
//          Eigen::Map<Eigen::MatrixXd>gradGear(grad.data()+offset, 2, 1);
//          offset += 2;
//          gradGear.setZero();
//          Gear_container.push_back(Gear);
//          gradGear_container.push_back(gradGear);
//      }
// //        //
//      Eigen::Map<const Eigen::VectorXd> Angles(x.data()+offset, opt->trajnum-1);
//      Eigen::Map<Eigen::VectorXd>gradAngles(grad.data()+offset, opt->trajnum-1);
//      gradAngles.setZero();
//
//        /*for relax end v?*/
//        // offset += opt->trajnum-1;
//        // Eigen::Map<const Eigen::MatrixXd> endV (x.data()+offset, 2, 1);
//        // Eigen::Map<Eigen::MatrixXd> gradendV (grad.data()+offset, 2, 1);
//        // gradendV.setZero();
//
//
//
//
//
//
//
     // Eigen::VectorXd gradt;
     Eigen::VectorXd obs_surround_feas_qvar_costs(6);
     obs_surround_feas_qvar_costs.setZero();

     for(int trajid = 0; trajid < opt->trajnum; trajid++){
         double smoo_cost;
        // 0 1 2 3:l_v 4:phi_dot 5:phi 

        //  Eigen::MatrixXd IniS,FinS;
        //  IniS = opt->iniState_container[trajid];
        //  FinS = opt->finState_container[trajid];

        //  if(trajid > 0){
        //      double theta = Angles[trajid-1];
        //      IniS.col(0) = Gear_container[trajid-1];
        //      IniS.col(1) = Eigen::Vector2d(-opt->non_sinv*cos(theta), -opt->non_sinv*sin(theta));
        //  }
        //  if(trajid < opt->trajnum-1){
        //      double theta = Angles[trajid];
        //      FinS.col(0) = Gear_container[trajid];
        //      FinS.col(1) = Eigen::Vector2d(opt->non_sinv*cos(theta), opt->non_sinv*sin(theta));
        //  }
// //            //relax end v?
// //            // if(trajid == opt->trajnum - 1){
// //            //   FinS.col(1) = endV;
// //            // }
// //
        opt->jerkOpt_S3NU_container[trajid].generate(P_container[trajid], T);
        // opt->jerkOpt_S3NU_container[trajid].initSmGradCost_zeros();
        opt->jerkOpt_S3NU_container[trajid].initSmGradCost();
        smoo_cost = opt->jerkOpt_S3NU_container[trajid].getTrajJerkCost();


       

        //  opt->jerkOpt_container[trajid].generate(P_container[trajid],T[trajid] / opt->piece_num_container[trajid],IniS,FinS);
        //  opt->jerkOpt_container[trajid].initSmGradCost(); // Smoothness cost
        //  smoo_cost = opt->jerkOpt_container[trajid].getTrajJerkCost();
        // std::cout<<" enter addPVAGradCost2CT_UN"<<std::endl;

         opt->addPVAGradCost2CT_UN( obs_surround_feas_qvar_costs, trajid, 0.0); // Time int cost
        //  std::cout<<" out addPVAGradCost2CT_UN"<<std::endl;

         //Get gradT gradC
         total_smcost += smoo_cost;
         penalty_cost +=  obs_surround_feas_qvar_costs.sum();
// //            // std::cout<<"Trajid: "<<trajid<<" penalty: "<<obs_surround_feas_qvar_costs.transpose()<<std::endl;
     }
// //
// //
     for(int trajid = 0; trajid < opt->trajnum; trajid++){
         double time_cost = 0.0;
        //  Eigen::Matrix<double,2,3> gradIni, gradFin;
        // std::cout<<" 1"<<std::endl;

         opt->jerkOpt_S3NU_container[trajid].calGrads_PT(); // gdt gdp gdhead gdtail
         //waypoint
        //  std::cout<<" 2"<<std::endl;

         gradP_container[trajid] = opt->jerkOpt_S3NU_container[trajid].get_gdP();
        //  std::cout<<" 3"<<std::endl;

         //init Fin
        //  gradIni = opt->jerkOpt_container[trajid].get_gdHead();
        //  gradFin = opt->jerkOpt_container[trajid].get_gdTail();
        //  if(opt->GearOpt){
        //      if(trajid > 0){
        //          double theta = Angles[trajid-1];
        //          gradGear_container[trajid-1] += gradIni.col(0);
        //          gradAngles[trajid-1] += gradIni.col(1).transpose() * Eigen::Vector2d(opt->non_sinv * sin(theta), -opt->non_sinv*cos(theta));
        //      }
        //      if(trajid < opt->trajnum-1){
        //          double theta = Angles[trajid];
        //          gradGear_container[trajid] += gradFin.col(0);
        //          gradAngles[trajid] += gradFin.col(1).transpose() * Eigen::Vector2d(-opt->non_sinv * sin(theta), opt->non_sinv*cos(theta));
        //      }
        //  }

         // if(trajid ==  opt->trajnum-1){
         //   gradendV += gradFin.col(1);
         // }
         Eigen::VectorXd gradt_vec = gradt;  // 将 gradt 映射转换为一个临时的 Eigen::VectorXd

         opt->VirtualTGradCost_UN(T,t,opt->jerkOpt_S3NU_container[trajid].get_gdT(),gradt_vec,time_cost);
         gradt = gradt_vec;  // 将修改后的值写回到 gradt
         // gradt[trajid] = 0.0;

         total_timecost += time_cost;
         // std::cout<<"gradp: \n"<<gradP_container[trajid] <<std::endl;
         // std::cout<<"gradeigen p: \n"<<opt->jerkOpt_container[trajid].get_gdP()<<std::endl;
     }
     
  //    std::cout<<"opt->wei_ref_: " <<opt->wei_ref_<<std::endl;
      for(int i = 0; i<opt->piece_num_container[0]-1; i++){
          // totoal_refcost += x_l[i]*x_l[i];
          grad_x_l[i]=-sin(opt->initInnerPts_ref_container[0](3,i))*gradP_container[0](0,i)+
          cos(opt->initInnerPts_ref_container[0](3,i))*gradP_container[0](1,i);
          // grad_x_l[i] += opt->wei_ref_*2*x_l[i];
      }
      grad_x_l.tail(piece_num_)=gradt;

      // totoal_refcost = totoal_refcost *opt->wei_ref_;
      // std::cout<<"totoal_refcost: " << totoal_refcost <<std::endl;
     
     // std::cout << "grad angle : "<< gradAngles.transpose() << "\n";
     // std::cout << "angle : "<< Angles.transpose() << "\n";
     // std::cout <<"endV: " <<endV.transpose()<<"\n";
     opt->RecordIterData(P_container);

      opt->recode_data_cost(total_smcost,obs_surround_feas_qvar_costs,total_timecost);
      opt->iter_num_ += 1;
     return total_smcost + total_timecost + penalty_cost 
     ;
}

void PolyTrajOptimizer::addPVAGradCost2CT_UN(Eigen::VectorXd &costs,  const int trajid, const double trajtime)
  {

    // output gradT gradC 
    int N = piece_num_container[trajid];
    Eigen::Vector2d outerNormal;
    Eigen::Vector2d sigma, dsigma, ddsigma, dddsigma, ddddsigma;
    double vel2_reci,vel2_reci_e,vel3_2_reci_e,acc2, cur2, cur, phi_dot;
    double latacc2;
    Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double step, alpha;
    Eigen::Matrix<double, 6, 2> gradPhiDotc,gradRefAc, gradRefVc,gradRefSc, gradViolaPc, gradViolaVc, gradViolaAc,gradViolaLatAc, gradViolaKc,gradViolaKLc,gradViolaKRc,gradViolaPhidotLc,gradViolaPhidotRc;
    double gradPhiDott, gradRefSt,gradRefAt,gradRefVt,gradViolaPt, gradViolaVt, gradViolaAt,gradViolaLatAt, gradViolaKt,gradViolaKLt,gradViolaKRt,gradViolaPhidotLt,gradViolaPhidotRt;
    double violaPos, violaVel, violaAcc, violaLatAcc, violaCur, violaCurL, violaCurR, violaDynamicObs, violaPhidotL, violaPhidotR;
    double violaPhidot, refSPenaD,refAccPenaD ,phiDotPenaD, refVelPenaD, violaPosPenaD, violaVelPenaD, violaAccPenaD, violaLatAccPenaD, violaCurPenaD, violaCurPenaDL, violaCurPenaDR,violaDynamicObsPenaD, violaPhidotPenaDL, violaPhidotPenaDR;
    double phiDotPena,refAccPena,refVelPena,refSPena,violaPosPena, violaVelPena, violaAccPena, violaLatAccPena, violaCurPena, violaCurPenaL, violaCurPenaR,violaDynamicObsPena, violaPhidotPenaL, violaPhidotPenaR;
    double phidot_denominator, phidot_nominator;

    double approxcur2, approxviolaCur,approxviolaCurPenaD,approxviolaCurPena;
    Eigen::Matrix<double, 6, 2> gradapproxViolaKc;
    double gradapproxViolaKt;
    
    std::vector<Eigen::MatrixXd> cfgHs = cfgHs_container[trajid];
    int singul_ = singul_container[trajid];
    double max_vel,max_cur,max_acc;
    if(singul_ > 0){
      max_vel = max_forward_vel;
      max_cur = max_forward_cur;
      max_acc = max_forward_acc;
    }
    else{
      max_vel = max_backward_vel;
      max_cur = max_backward_cur;
      max_acc = max_backward_acc;
    }


    double omg;
    int i_dp = 0; // the index of constrain points
    costs.setZero();
    double z_h0, z_h1, z_h2, z_h3, z_h4;
    double z1, z2, z3;
    double n1, n2, n3, n4, n5, n6;
    Eigen::Matrix2d ego_R, help_R;
    /*debug*/
    cos_points.clear();
    debug_hPolys.clear();
    key_points.clear();
    std::vector<int> cosindex;
    // int innerLoop;
    double t = 0;

    //debug
    double velcost=0.0;
    double acccost=0.0;
    double latacccost = 0.0;
    double curcost=0.0;
    double phidotcost = 0.0;
    int pointid = -1;




    for (int i = 0; i < N; ++i)
    {
      int K;
      if(i==0 || i==N-1){
        K = destraj_resolution_;
      }
      else{
        K = traj_resolution_;
      }
      const Eigen::Matrix<double, 6, 2> &c = jerkOpt_S3NU_container[trajid].getCoeffs().block<6, 2>(i * 6, 0);
      const double piece_time = jerkOpt_S3NU_container[trajid].getT()(i);
      step = piece_time / K; // T_i /k
      s1 = 0.0;
      // innerLoop = K;
      for (int j = 0; j <= K; ++j)
      {
        s2 = s1 * s1;
        s3 = s2 * s1;
        s4 = s2 * s2;
        s5 = s4 * s1;
        beta0 << 1.0, s1, s2, s3, s4, s5;
        beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
        beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
        beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
        beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120 * s1;
        alpha = 1.0 / K * j;
        
        //update s1 for the next iteration
        s1 += step;
        pointid++;

        sigma = c.transpose() * beta0;
        dsigma = c.transpose() * beta1;
        ddsigma = c.transpose() * beta2;
        dddsigma = c.transpose() * beta3;
        ddddsigma = c.transpose() * beta4;
         
        // ctrl_points_.col(i_dp) = sigma;
        omg = (j == 0 || j == K) ? 0.5 : 1.0;

        // some help values
        
        z_h0 = dsigma.norm();
        z_h1 = ddsigma.transpose() * dsigma;
        z_h2 = dddsigma.transpose() * dsigma;
        z_h3 = ddsigma.transpose() * B_h * dsigma;

        n1 = z_h0;
        n2 = n1 * n1;
        n3 = n2 * n1;
        n4 = n2 * n2;
        n5 = n3 * n2;
        n6 = n3 * n3;
      
        z1 = dddsigma.transpose() * B_h * dsigma;
        z2 = ddsigma.transpose() * B_h * dsigma;
        z3 = dsigma.transpose() * ddsigma;  
        
        if (j != K || (j == K && i == N - 1))
        {
          ++i_dp;
        }


        // add cost z_h0 = ||v||
        if ( z_h0 < 1e-4 || (j==0&&i==0) || (i==N-1&&j==K))
        {
          continue;
        }
        //avoid siguality

        vel2_reci = 1.0 / (z_h0 * z_h0);
        vel2_reci_e = 1.0 / (z_h0 * z_h0+epis);
        vel3_2_reci_e = vel2_reci_e * sqrt(vel2_reci_e);
        z_h0 = 1.0 / z_h0;

        z_h4 = z_h1 * vel2_reci;
        violaVel = 1.0 / vel2_reci - max_vel * max_vel;
        acc2 = z_h1 * z_h1 * vel2_reci;
        latacc2 = z_h3 * z_h3 * vel2_reci;
        cur2 = z_h3 * z_h3 * (vel2_reci_e * vel2_reci_e * vel2_reci_e);
        cur = z_h3 * vel3_2_reci_e;
        violaAcc = acc2 - max_acc * max_acc;
        violaLatAcc = latacc2 - max_latacc_ * max_latacc_;

        phidot_denominator = n6 + L_ * L_ * z2 * z2;
        phidot_nominator = L_ * (n3 * z1 - 3 * z2 * z3 * n1);
        phi_dot = phidot_nominator / phidot_denominator; // S/M

        //@hzc: add feasibility with curvature
        violaCur = cur2 - max_cur * max_cur;
        violaCurL = cur-max_cur;
        violaCurR = -cur-max_cur;


        ego_R << dsigma(0), -dsigma(1),
                 dsigma(1), dsigma(0);
        ego_R = singul_ * ego_R * z_h0;

        Eigen::Matrix2d temp_a, temp_v;
        temp_a << ddsigma(0), -ddsigma(1),
                  ddsigma(1), ddsigma(0);
        temp_v << dsigma(0), -dsigma(1),
                  dsigma(1), dsigma(0);
        Eigen::Matrix2d R_dot = singul_ * (temp_a * z_h0 - temp_v * vel2_reci * z_h0 * z_h1);

        for(auto le : vec_le_)
        {
          Eigen::Vector2d bpt = sigma + ego_R * le;

          Eigen::Matrix2d temp_l_Bl;
          temp_l_Bl << le(0), -le(1),
                       le(1), le(0);          

          int corr_k = cfgHs[pointid].cols();

          for(int k = 0; k < corr_k; k++)
          {
            outerNormal = cfgHs[pointid].col(k).head<2>();
            violaPos = outerNormal.dot(bpt - cfgHs[pointid].col(k).tail<2>());

            if(violaPos > 0)
            {
              positiveSmoothedL1(violaPos, violaPosPena, violaPosPenaD);
              
              gradViolaPc = beta0 * outerNormal.transpose() + 
                            beta1 * outerNormal.transpose() * (singul_ * temp_l_Bl * z_h0 - ego_R * le * dsigma.transpose() * vel2_reci);
              
              gradViolaPt = alpha * outerNormal.transpose() * (dsigma + R_dot * le);

              jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_obs_ * violaPosPenaD * gradViolaPc;
              jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_obs_ * (violaPosPenaD * gradViolaPt * step + violaPosPena / K);

              costs(0) += omg * step * wei_obs_ * violaPosPena; // cost is the same
            }
          }
        }
        
        // ---------------------surrounding vehicle avoidance
        double gradt, grad_prev_t, costp;
        Eigen::Vector2d gradp, gradp2;


        // if(surroundGradCostP(i_dp, t + step * j, sigma, dsigma, gradp, gradt, grad_prev_t, costp))
        // {

        //signed dist 
        // use it! @hzc


        if(surround_trajs_!=NULL){
          costs(1) += dynamicObsGradCostP(omg,step,t + step * j,beta0,beta1,alpha,i,K,sigma,dsigma,ddsigma,ego_R,help_R,trajid,trajtime);
        }

//       //min(phi_dot)
        if(wei_u_phi_dot_  > 0.001){
          violaPhidotR = fabs(phi_dot);
          double phidot_dection=0.0;
          if(phi_dot >= 0){
            phidot_dection=1.0;
          }else{
            phidot_dection=-1.0;
          }
          positiveSmoothedL1(violaPhidotR, violaPhidotPenaR, violaPhidotPenaDR);
          Eigen::Vector2d partial_S_over_partial_dsigma
              = L_ * (n3 * B_h.transpose() * dddsigma + 3 * z1 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z3 * n1 - 3 * z2 * ddsigma * n1 - 3 * z2 * z3 * dsigma / n1);
          Eigen::Vector2d partial_M_over_partial_dsigma 
              = 6 * n4 * dsigma + 2 * L_ * L_ * z2 * B_h.transpose() * ddsigma;
          Eigen::Vector2d partial_S_over_partial_ddsigma
              = -3 * L_ * n1 * (B_h * dsigma * z3 + z2 * dsigma);
          Eigen::Vector2d partial_M_over_partial_ddsigma
              = 2 * L_ * L_ * z2 * B_h * dsigma;
          Eigen::Vector2d partial_S_over_partial_dddsigma
              = L_ * n3 * B_h * dsigma;

          Eigen::Vector2d partial_phi_dot_over_partial_dsigma
              = (partial_S_over_partial_dsigma * phidot_denominator - partial_M_over_partial_dsigma * phidot_nominator) / pow(phidot_denominator, 2);
          Eigen::Vector2d partial_phi_dot_over_partial_ddsigma
              = (partial_S_over_partial_ddsigma * phidot_denominator - partial_M_over_partial_ddsigma * phidot_nominator) / pow(phidot_denominator, 2);
          Eigen::Vector2d partial_phi_dot_over_partial_dddsigma
              = partial_S_over_partial_dddsigma / phidot_denominator;

          gradViolaPhidotRc = phidot_dection* beta1 * partial_phi_dot_over_partial_dsigma.transpose()
                            +phidot_dection* beta2 * partial_phi_dot_over_partial_ddsigma.transpose()
                            +phidot_dection* beta3 * partial_phi_dot_over_partial_dddsigma.transpose();
          gradViolaPhidotRt = phidot_dection* alpha * (partial_phi_dot_over_partial_dsigma.transpose() * ddsigma
                                      +partial_phi_dot_over_partial_ddsigma.transpose() * dddsigma
                                      +partial_phi_dot_over_partial_dddsigma.transpose() * ddddsigma)(0, 0);

          jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_u_phi_dot_ * 1.0 * violaPhidotPenaDR * gradViolaPhidotRc;
          jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_u_phi_dot_ * 1.0 * (violaPhidotPenaDR * gradViolaPhidotRt * step + violaPhidotPenaR / K);
          costs(4) += omg * step * wei_u_phi_dot_ * 1.0 * violaPhidotPenaR;
          phidotcost+=omg * step * wei_u_phi_dot_ * 1.0 * violaPhidotPenaR;
        }




        if( i < N - 1 && i < N - 1){

          //error between innerpt and reference line
            if(wei_ref_ > 0.001 ){
            Eigen::VectorXd p_e = initInnerPts_ref_container[0].col(i);
            Eigen::Vector2d p_r(p_e(1),p_e(2));

            Eigen::Vector2d N_r(-sin(p_e(3)),cos(p_e(3)));
            
            double e_s = ((sigma - p_r).transpose()*N_r)(0,0);
            Eigen::Vector2d partial_e_v_over_partial_sigma 
            = N_r;
            refSPena = e_s * e_s;
            refSPenaD = 2 * e_s;
            gradRefSc = beta0 * partial_e_v_over_partial_sigma.transpose();
            gradRefSt = (partial_e_v_over_partial_sigma.transpose()*dsigma)(0,0);
            jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += wei_ref_ * piece_time  * refSPenaD * gradRefSc;
            jerkOpt_S3NU_container[trajid].get_gdT()(i) += wei_ref_ * (piece_time * refSPenaD * gradRefSt + refSPena);
            costs(3) += wei_ref_* refSPena;
            }

          //error v between innerpt and reference line
          if(wei_ref_v_ > 0.001){

            // de/dt
                  // Eigen::VectorXd p_e = initInnerPts_ref_container[0].col(i);
                  // Eigen::Vector2d N_r(-sin(p_e(3)),cos(p_e(3)));
                  // double e_v = dsigma.transpose()*N_r;
                  // Eigen::Vector2d partial_e_v_over_partial_dsigma 
                  //                 = N_r;
                  // refVelPena = e_v * e_v;
                  // refVelPenaD = 2 * e_v;
                  // gradRefVc = beta1 * partial_e_v_over_partial_dsigma.transpose();
                  // gradRefVt = (partial_e_v_over_partial_dsigma.transpose() * ddsigma)(0,0);
                  //   jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += wei_ref_v_ * piece_time * refVelPenaD * gradRefVc;
                  //   jerkOpt_S3NU_container[trajid].get_gdT()(i) += wei_ref_v_ * (piece_time * refVelPenaD * gradRefVt + refVelPena);
                  // costs(3) += wei_ref_v_* refVelPena;

            // de/ds

                  Eigen::VectorXd p_e = initInnerPts_ref_container[0].col(i);
                  Eigen::Vector2d N_r(-sin(p_e(3)),cos(p_e(3)));
                  Eigen::Vector2d T_r( cos(p_e(3)),sin(p_e(3)));
                  Eigen::Vector2d p_r(p_e(1),p_e(2));
                  double e = (sigma - p_r).norm();
                  double k_r = p_e(4);


                  double e_v = (1-k_r*e)*(dsigma.transpose()*N_r)(0,0)/(dsigma.transpose()*T_r)(0,0);

                  // double e_v = dsigma.transpose()*N_r;
                  Eigen::Vector2d partial_e_v_over_partial_sigma 
                                  = - k_r*dsigma.transpose()*N_r*N_r/(dsigma.transpose()*T_r)(0,0);
                    Eigen::Vector2d partial_e_v_over_partial_dsigma 
                            = (1 - k_r * e) * (
                              N_r * (dsigma.transpose() * T_r)(0, 0) 
                              - T_r * (dsigma.transpose() * N_r)(0, 0)) 
                              / pow((dsigma.transpose() * T_r)(0, 0), 2);
                  
                  refVelPena = e_v * e_v;
                  refVelPenaD = 2 * e_v;

                  gradRefVc = beta0 * partial_e_v_over_partial_sigma.transpose()
                              + beta1 * partial_e_v_over_partial_dsigma.transpose();

                  gradRefVt = (partial_e_v_over_partial_sigma.transpose() * dsigma
                             + partial_e_v_over_partial_dsigma.transpose()*ddsigma)(0,0);

                    jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += wei_ref_v_ * piece_time * refVelPenaD * gradRefVc;
                    jerkOpt_S3NU_container[trajid].get_gdT()(i) += wei_ref_v_ * (piece_time * refVelPenaD * gradRefVt + refVelPena);
                    costs(3) += wei_ref_v_* refVelPena;
          }


          // //error a between innerpt and reference line
          //     if(wei_ref_a_ > 0.001){

          //             Eigen::VectorXd p_e = initInnerPts_ref_container[0].col(i);
          //             Eigen::Vector2d N_r(-sin(p_e(3)),cos(p_e(3)));
      
          //             Eigen::Vector2d p_r(p_e(1),p_e(2));
      
          //             double e = (sigma - p_r).norm();
          //             double k_r = p_e(4);
      
          //             double e_a = ddsigma.transpose()*N_r-(dsigma[0] * dsigma[0] + dsigma[1] * dsigma[1])*k_r/(1-k_r*e);
      
          //             Eigen::Vector2d partial_e_a_over_partial_sigma 
          //                               = -k_r * k_r * (dsigma.transpose() * dsigma)(0,0) * (sigma - p_r) / (e * pow(1 - k_r * e, 2));
          //             Eigen::Vector2d partial_e_a_over_partial_dsigma 
          //                               = - dsigma * (2*k_r)/(1-k_r*e);
          //             Eigen::Vector2d partial_e_a_over_partial_ddsigma 
          //                               = N_r;
          //             refAccPena = e_a * e_a;
          //             refAccPenaD = 2 * e_a;
          //             gradRefAc = beta0 * partial_e_a_over_partial_sigma.transpose()+
          //                         beta1 * partial_e_a_over_partial_dsigma.transpose() +
          //                         beta2 * partial_e_a_over_partial_ddsigma.transpose();
          //             gradRefAt = (partial_e_a_over_partial_sigma.transpose()*dsigma +
          //                         partial_e_a_over_partial_dsigma.transpose()*ddsigma +
          //                         partial_e_a_over_partial_ddsigma.transpose()*dddsigma)(0,0);
          //             jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += wei_ref_a_ * piece_time  * refAccPenaD * gradRefAc;
          //             jerkOpt_S3NU_container[trajid].get_gdT()(i) += wei_ref_a_ * (piece_time * refAccPenaD * gradRefAt + refAccPena);
          //             costs(3) += wei_ref_a_* refVelPena;
          //   }

        }



          if (violaVel > 0.0)
        {
          positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);

          gradViolaVc = 2.0 * beta1 * dsigma.transpose(); // 6*2
          gradViolaVt = 2.0 * alpha * z_h1;               // 1*1
          jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaVelPenaD * gradViolaVc;
          jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_feas_ * (violaVelPenaD * gradViolaVt * step + violaVelPena / (K));
          costs(2) += omg * step * wei_feas_ * violaVelPena;
          velcost+=omg * step * wei_feas_ * violaVelPena;

        }
        
        if (violaAcc > 0.0)
        {
          positiveSmoothedL1(violaAcc, violaAccPena, violaAccPenaD);
          gradViolaAc = 2.0 * beta1 * (z_h4 * ddsigma.transpose() - z_h4 * z_h4 * dsigma.transpose()) +
                        2.0 * beta2 * z_h4 * dsigma.transpose(); // 6*2
          gradViolaAt = 2.0 * alpha * (z_h4 * (ddsigma.squaredNorm() + z_h2) - z_h4 * z_h4 * z_h1);
          jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaAccPenaD * gradViolaAc;
          jerkOpt_S3NU_container[trajid].get_gdT()(i)  += omg * wei_feas_ * (violaAccPenaD * gradViolaAt * step + violaAccPena /  (K));
          costs(2) += omg * step * wei_feas_ * violaAccPena;
          acccost += omg * step * wei_feas_ * violaAccPena;
        } 
        // if (violaLatAcc > 0.0)
        // {
        //   positiveSmoothedL1(violaLatAcc, violaLatAccPena, violaLatAccPenaD);
        //   gradViolaLatAc = 2.0 * beta1 * (z_h3 * vel2_reci * ddsigma.transpose() * B_h - z_h3 * vel2_reci  * z_h3 * vel2_reci  * dsigma.transpose()) +
        //                 2.0 * beta2 * z_h3  * vel2_reci * dsigma.transpose() * B_h.transpose(); // 6*2
        //   gradViolaLatAt = 2.0 * alpha * (z_h3  * vel2_reci * z1
        //                     -z_h3 * vel2_reci  * z_h3 * vel2_reci * dddsigma.transpose()*B_h*dsigma);

        //   jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaLatAccPenaD * gradViolaLatAc;
        //   jerkOpt_container[trajid].get_gdT()  += omg * wei_feas_ * (violaLatAccPenaD * gradViolaLatAt * step + violaLatAccPena / K);
        //   costs(2) += omg * step * wei_feas_ * violaLatAccPena;
        //   latacccost += omg * step * wei_feas_ * violaLatAccPena;
        // } 

       
        /*violaCurL = cur-max_cur_;
        violaCurR = -cur-max_cur_;*/

        if(violaCurL > 0.0){
          positiveSmoothedL1(violaCurL, violaCurPenaL, violaCurPenaDL);
          //@hzc
          gradViolaKLc = beta1 * (vel3_2_reci_e * ddsigma.transpose()*B_h - 3 * vel3_2_reci_e * vel2_reci_e * z_h3 * dsigma.transpose()) 
                         + beta2 * vel3_2_reci_e * dsigma.transpose() * B_h.transpose(); // 6*2
          gradViolaKLt  = alpha*vel3_2_reci_e*(dddsigma.transpose()*B_h*dsigma-3*vel2_reci_e*z_h3*z_h1);
          jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * 10.0 * violaCurPenaDL * gradViolaKLc;
          jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_feas_ * 10.0 * (violaCurPenaDL * gradViolaKLt * step + violaCurPenaL /  (K));
          costs(2) += omg * step * wei_feas_ * 10.0 * violaCurPenaL;
          curcost+=omg * step * wei_feas_ * 10.0 * violaCurPenaL;
        }
        if(violaCurR > 0.0){
          positiveSmoothedL1(violaCurR, violaCurPenaR, violaCurPenaDR);
          //@hzc
          gradViolaKRc = -(beta1 * (vel3_2_reci_e * ddsigma.transpose()*B_h - 3 * vel3_2_reci_e * vel2_reci_e * z_h3 * dsigma.transpose()) 
                         + beta2 * vel3_2_reci_e * dsigma.transpose() * B_h.transpose()); // 6*2
          gradViolaKRt  = -(alpha*vel3_2_reci_e*(dddsigma.transpose()*B_h*dsigma-3*vel2_reci_e*z_h3*z_h1));
          jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * 10.0 * violaCurPenaDR * gradViolaKRc;
          jerkOpt_S3NU_container[trajid].get_gdT()(i)  += omg * wei_feas_ * 10.0 * (violaCurPenaDR * gradViolaKRt * step + violaCurPenaR /  (K));
          costs(2) += omg * step * wei_feas_ * 10.0 * violaCurPenaR;
          curcost+=omg * step * wei_feas_ * 10.0 * violaCurPenaR;
        }



        // if(violaCurL > 0.0){
        //   positiveSmoothedL1(violaCurL, violaCurPenaL, violaCurPenaDL);
        //   //@hzc
        //   gradViolaKLc = beta1 * (vel3_2_reci_e * ddsigma.transpose()*B_h - 3 * vel3_2_reci_e * vel2_reci_e * z_h3 * dsigma.transpose()) 
        //                  + beta2 * vel3_2_reci_e * dsigma.transpose() * B_h.transpose(); // 6*2
        //   gradViolaKLt  = alpha*vel3_2_reci_e*(dddsigma.transpose()*B_h*dsigma-3*vel2_reci_e*z_h3*z_h1);
        //   jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_ref_a * 10.0 * violaCurPenaDL * gradViolaKLc;
        //   jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_ref_a * 10.0 * (violaCurPenaDL * gradViolaKLt * step + violaCurPenaL /  (K));
        //   costs(2) += omg * step * wei_ref_a * 10.0 * violaCurPenaL;
        //   curcost+=omg * step * wei_ref_a * 10.0 * violaCurPenaL;
        // }
        // if(violaCurR > 0.0){
        //   positiveSmoothedL1(violaCurR, violaCurPenaR, violaCurPenaDR);
        //   //@hzc
        //   gradViolaKRc = -(beta1 * (vel3_2_reci_e * ddsigma.transpose()*B_h - 3 * vel3_2_reci_e * vel2_reci_e * z_h3 * dsigma.transpose()) 
        //                  + beta2 * vel3_2_reci_e * dsigma.transpose() * B_h.transpose()); // 6*2
        //   gradViolaKRt  = -(alpha*vel3_2_reci_e*(dddsigma.transpose()*B_h*dsigma-3*vel2_reci_e*z_h3*z_h1));
        //   jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_ref_a * 10.0 * violaCurPenaDR * gradViolaKRc;
        //   jerkOpt_S3NU_container[trajid].get_gdT()(i)  += omg * wei_ref_a * 10.0 * (violaCurPenaDR * gradViolaKRt * step + violaCurPenaR /  (K));
        //   costs(2) += omg * step * wei_ref_a * 10.0 * violaCurPenaR;
        //   curcost+=omg * step * wei_ref_a * 10.0 * violaCurPenaR;
        // }



        // violaPhidotL = phi_dot - max_phidot_;
        // violaPhidotR = -phi_dot - max_phidot_;




        // phiDotPena= phi_dot * phi_dot;
        // phiDotPenaD = 2.0 * phi_dot;
        // Eigen::Vector2d partial_S_over_partial_dsigma
        //       = L_ * (n3 * B_h.transpose() * dddsigma + 3 * z1 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z3 * n1 - 3 * z2 * ddsigma - 3 * z2 * z3 * dsigma / z1);
        //   Eigen::Vector2d partial_M_over_partial_dsigma 
        //       = 6 * n4 * dsigma + 2 * L_ * L_ * z2 * B_h.transpose() * ddsigma;
        //   Eigen::Vector2d partial_S_over_partial_ddsigma
        //       = -3 * L_ * n1 * (B_h * dsigma * z3 + z2 * dsigma);
        //   Eigen::Vector2d partial_M_over_partial_ddsigma
        //       = 2 * L_ * L_ * z2 * B_h * dsigma;
        //   Eigen::Vector2d partial_S_over_partial_dddsigma
        //       = L_ * n3 * B_h * dsigma;
        //   Eigen::Vector2d partial_phi_dot_over_partial_dsigma
        //       = (partial_S_over_partial_dsigma * phidot_denominator - partial_M_over_partial_dsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_ddsigma
        //       = (partial_S_over_partial_ddsigma * phidot_denominator - partial_M_over_partial_ddsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_dddsigma
        //       = partial_S_over_partial_dddsigma / phidot_denominator;

        // gradPhiDotc = beta1 * partial_phi_dot_over_partial_dsigma.transpose()
        //                       + beta2 * partial_phi_dot_over_partial_ddsigma.transpose()
        //                       + beta3 * partial_phi_dot_over_partial_dddsigma.transpose();
        // gradPhiDott= alpha * (partial_phi_dot_over_partial_dsigma.transpose() * ddsigma
        //                     +partial_phi_dot_over_partial_ddsigma.transpose() * dddsigma
        //                      +partial_phi_dot_over_partial_dddsigma.transpose() * ddddsigma)(0, 0);
        //                     //  double wei_phi_dot_=200.0;
        // jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step *  wei_ref_a* phiDotPenaD * gradPhiDotc;
        // jerkOpt_S3NU_container[trajid].get_gdT()(i)  += omg  * wei_ref_a * ( phiDotPenaD * gradPhiDott * step + phiDotPena / K);
        // costs(4) += omg * step * wei_ref_a * phiDotPena;

        // if(violaPhidotL > 0.0)
        // {
        //   positiveSmoothedL1(violaPhidotL, violaPhidotPenaL, violaPhidotPenaDL);
        //   Eigen::Vector2d partial_S_over_partial_dsigma
        //       = L_ * (n3 * B_h.transpose() * dddsigma + 3 * z1 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z3 * n1 - 3 * z2 * ddsigma - 3 * z2 * z3 * dsigma / z1);
        //   Eigen::Vector2d partial_M_over_partial_dsigma 
        //       = 6 * n4 * dsigma + 2 * L_ * L_ * z2 * B_h.transpose() * ddsigma;
        //   Eigen::Vector2d partial_S_over_partial_ddsigma
        //       = -3 * L_ * n1 * (B_h * dsigma * z3 + z2 * dsigma);
        //   Eigen::Vector2d partial_M_over_partial_ddsigma
        //       = 2 * L_ * L_ * z2 * B_h * dsigma;
        //   Eigen::Vector2d partial_S_over_partial_dddsigma
        //       = L_ * n3 * B_h * dsigma;

        //   Eigen::Vector2d partial_phi_dot_over_partial_dsigma
        //       = (partial_S_over_partial_dsigma * phidot_denominator - partial_M_over_partial_dsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_ddsigma
        //       = (partial_S_over_partial_ddsigma * phidot_denominator - partial_M_over_partial_ddsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_dddsigma
        //       = partial_S_over_partial_dddsigma / phidot_denominator;

        //   gradViolaPhidotLc = beta1 * partial_phi_dot_over_partial_dsigma.transpose()
        //                     + beta2 * partial_phi_dot_over_partial_ddsigma.transpose()
        //                     + beta3 * partial_phi_dot_over_partial_dddsigma.transpose();
        //   gradViolaPhidotLt = alpha * (partial_phi_dot_over_partial_dsigma.transpose() * ddsigma
        //                               +partial_phi_dot_over_partial_ddsigma.transpose() * dddsigma
        //                               +partial_phi_dot_over_partial_dddsigma.transpose() * ddddsigma)(0, 0);

        //   jerkOpt_S3NU_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_ref_a * 1.0 * violaPhidotPenaDL * gradViolaPhidotLc;
        //   jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_ref_a * 1.0 * (violaPhidotPenaDL * gradViolaPhidotLt * step + violaPhidotPenaL / K);
        //   costs(4) += omg * step * wei_ref_a * 1.0 * violaPhidotPenaL;
        //   phidotcost+=omg * step * wei_ref_a * 1.0 * violaPhidotPenaL;
        // }

        // if(violaPhidotR > 0.0)
        // {
        //   positiveSmoothedL1(violaPhidotR, violaPhidotPenaR, violaPhidotPenaDR);
        //   Eigen::Vector2d partial_S_over_partial_dsigma
        //       = L_ * (n3 * B_h.transpose() * dddsigma + 3 * z1 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z3 * n1 - 3 * z2 * ddsigma * n1 - 3 * z2 * z3 * dsigma / n1);
        //   Eigen::Vector2d partial_M_over_partial_dsigma 
        //       = 6 * n4 * dsigma + 2 * L_ * L_ * z2 * B_h.transpose() * ddsigma;
        //   Eigen::Vector2d partial_S_over_partial_ddsigma
        //       = -3 * L_ * n1 * (B_h * dsigma * z3 + z2 * dsigma);
        //   Eigen::Vector2d partial_M_over_partial_ddsigma
        //       = 2 * L_ * L_ * z2 * B_h * dsigma;
        //   Eigen::Vector2d partial_S_over_partial_dddsigma
        //       = L_ * n3 * B_h * dsigma;

        //   Eigen::Vector2d partial_phi_dot_over_partial_dsigma
        //       = (partial_S_over_partial_dsigma * phidot_denominator - partial_M_over_partial_dsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_ddsigma
        //       = (partial_S_over_partial_ddsigma * phidot_denominator - partial_M_over_partial_ddsigma * phidot_nominator) / pow(phidot_denominator, 2);
        //   Eigen::Vector2d partial_phi_dot_over_partial_dddsigma
        //       = partial_S_over_partial_dddsigma / phidot_denominator;

        //   gradViolaPhidotRc = -beta1 * partial_phi_dot_over_partial_dsigma.transpose()
        //                     - beta2 * partial_phi_dot_over_partial_ddsigma.transpose()
        //                     - beta3 * partial_phi_dot_over_partial_dddsigma.transpose();
        //   gradViolaPhidotRt = -alpha * (partial_phi_dot_over_partial_dsigma.transpose() * ddsigma
        //                               +partial_phi_dot_over_partial_ddsigma.transpose() * dddsigma
        //                               +partial_phi_dot_over_partial_dddsigma.transpose() * ddddsigma)(0, 0);

        //   // jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_ref_a * 1.0 * violaPhidotPenaDR * gradViolaPhidotRc;
        //   jerkOpt_S3NU_container[trajid].get_gdT()(i) += omg * wei_ref_a * 1.0 * (violaPhidotPenaDR * gradViolaPhidotRt * step + violaPhidotPenaR / K);
        //   costs(4) += omg * step * wei_ref_a * 1.0 * violaPhidotPenaR;
        //   phidotcost+=omg * step * wei_ref_a * 1.0 * violaPhidotPenaR;
        // }
      }
      t += jerkOpt_container[trajid].getDt();
    }


  }

  void PolyTrajOptimizer::RecordIterData(const std::vector<Eigen::Map<const Eigen::MatrixXd>> &P_container){
    int piece_num_ = piece_num_container[0];
    if(!is_record_data_){
      return;
    }
      if(iter_num_ % 50 == 0){
        std::cout << "Recording data for iteration: " << iter_num_ << std::endl;
        std::string traj_file = workspace_ + "src/npct_Planner/data_plot/data_raw/data_iter/npct_iter_" + std::to_string(iter_num_) + ".csv";
        std::cout << "write dir: " << traj_file << std::endl;
        std::remove(traj_file.c_str());
        std::ofstream csvfile_out;
        csvfile_out.open(traj_file,
          std::ios::out | std::ios::app );
          if(csvfile_out.is_open()) {
            uint16_t row_count = 0;

            plan_utils::Trajectory traj = jerkOpt_S3NU_container[0].getTraj(1);
            double total_duration = traj.getTotalDuration();
            for (double t = 0;t <= total_duration ; t += 0.1) {
              Eigen::Vector2d inner_point = traj.getPos(t);
              double npct_t = t;
              double npct_x = inner_point(0);
              double npct_y = inner_point(1);
              double npct_v = traj.getVel(t);
              double npct_a = traj.getAcc(t);
              double heading = traj.getAngle(t);
              double p1_x, p1_y, p2_x, p2_y, p3_x, p3_y, p4_x, p4_y;
              double length = 4.88;
              double width = 1.90;
              double d_cr = 1.015;
              Eigen::Matrix2d R;
              R << cos(heading),-sin(heading),
                    sin(heading),cos(heading);

                    // Eigen::Vector2d rear_center = inner_point;
                    // Eigen::Vector2d front_center = inner_point + R * Eigen::Vector2d(length, 0);
             Eigen::Vector2d rear_left = inner_point + R * Eigen::Vector2d(length/2.0+d_cr, width / 2.0);
             Eigen::Vector2d rear_right = inner_point + R * Eigen::Vector2d(length/2.0+d_cr, -width / 2.0);
             Eigen::Vector2d front_left = inner_point + R * Eigen::Vector2d(-length/2.0+d_cr, -width / 2.0);
             Eigen::Vector2d front_right = inner_point + R * Eigen::Vector2d(-length/2.0+d_cr, width / 2.0);

             p1_x = rear_left(0);
             p1_y = rear_left(1);
            p2_x = rear_right(0);
            p2_y = rear_right(1);
            p3_x = front_left(0);
            p3_y = front_left(1);
            p4_x = front_right(0);
            p4_y = front_right(1);

              if(row_count < piece_num_){
                double piece_duration = traj.getDurations()(row_count);
                if(row_count < piece_num_-1){
                  double Inner_x = P_container[0](0,row_count);
                  double Inner_y = P_container[0](1,row_count);
                  double Inner_ref_x = initInnerPts_ref_container[0](1,row_count);
                  double Inner_ref_y = initInnerPts_ref_container[0](2,row_count);
                  if(row_count == 0 ){
                    csvfile_out << "index" << ","
                              << "piece_duration" << ","
                              << "inner_x" << ","
                              << "inner_y" << ","
                              << "inner_ref_x" << ","
                              << "inner_ref_y" << ","
                              << "npct_t" << ","
                              << "npct_x" << ","
                              << "npct_y" << ","
                              << "npct_v" << ","
                              << "npct_a" << ","
                              << "p1_x" << ","
                              << "p1_y" << ","
                              << "p2_x" << ","
                              << "p2_y" << ","
                              << "p3_x" << ","
                              << "p3_y" << ","
                              << "p4_x" << ","
                              << "p4_y" << ","
                              << std::endl;
                  }
                  csvfile_out << row_count << ","
                            << std::fixed << std::setprecision(6) << piece_duration << ","
                            << std::fixed << std::setprecision(6) << Inner_x << ","
                            << std::fixed << std::setprecision(6) << Inner_y << ","
                            << std::fixed << std::setprecision(6) << Inner_ref_x << ","
                            << std::fixed << std::setprecision(6) << Inner_ref_y << ","
                            << std::fixed << std::setprecision(6) << npct_t << ","
                            << std::fixed << std::setprecision(6) << npct_x<< ","
                            << std::fixed << std::setprecision(6) << npct_y<< ","
                            << std::fixed << std::setprecision(6) << npct_v<< ","
                            << std::fixed << std::setprecision(6) << npct_a<< ","
                            << std::fixed << std::setprecision(6) << p1_x<< ","
                            << std::fixed << std::setprecision(6) << p1_y<< ","
                            << std::fixed << std::setprecision(6) << p2_x<< ","
                            << std::fixed << std::setprecision(6) << p2_y<< ","
                            << std::fixed << std::setprecision(6) << p3_x<< ","
                            << std::fixed << std::setprecision(6) << p3_y<< ","
                            << std::fixed << std::setprecision(6) << p4_x<< ","
                            << std::fixed << std::setprecision(6) << p4_y<< ","<<std::endl;
                }else {
                  csvfile_out << row_count << ","
                  << std::fixed << std::setprecision(6) << piece_duration << ","
                  << std::fixed << std::setprecision(6) << "" << ","
                  << std::fixed << std::setprecision(6) << "" << ","
                  << std::fixed << std::setprecision(6) << "" << ","
                  << std::fixed << std::setprecision(6) << "" << ","
                  << std::fixed << std::setprecision(6) << npct_t << ","
                  << std::fixed << std::setprecision(6) << npct_x<< ","
                  << std::fixed << std::setprecision(6) << npct_y<< ","
                  << std::fixed << std::setprecision(6) << npct_v<< ","
                  << std::fixed << std::setprecision(6) << npct_a<< ","
                  << std::fixed << std::setprecision(6) << p1_x<< ","
                  << std::fixed << std::setprecision(6) << p1_y<< ","
                  << std::fixed << std::setprecision(6) << p2_x<< ","
                  << std::fixed << std::setprecision(6) << p2_y<< ","
                  << std::fixed << std::setprecision(6) << p3_x<< ","
                  << std::fixed << std::setprecision(6) << p3_y<< ","
                  << std::fixed << std::setprecision(6) << p4_x<< ","
                  << std::fixed << std::setprecision(6) << p4_y<< ","<<std::endl;
                }
              }else{
                csvfile_out << row_count << ","
                          << std::fixed << std::setprecision(6) << "" << ","
                          << std::fixed << std::setprecision(6) << "" << ","
                          << std::fixed << std::setprecision(6) << "" << ","
                          << std::fixed << std::setprecision(6) << "" << ","
                          << std::fixed << std::setprecision(6) << "" << ","
                          << std::fixed << std::setprecision(6) << npct_t << ","
                          << std::fixed << std::setprecision(6) << npct_x<< ","
                          << std::fixed << std::setprecision(6) << npct_y<< ","
                          << std::fixed << std::setprecision(6) << npct_v<< ","
                          << std::fixed << std::setprecision(6) << npct_a<< ","

                          << std::fixed << std::setprecision(6) << p1_x<< ","
                          << std::fixed << std::setprecision(6) << p1_y<< ","
                          << std::fixed << std::setprecision(6) << p2_x<< ","
                          << std::fixed << std::setprecision(6) << p2_y<< ","
                          << std::fixed << std::setprecision(6) << p3_x<< ","
                          << std::fixed << std::setprecision(6) << p3_y<< ","
                          << std::fixed << std::setprecision(6) << p4_x<< ","
                          << std::fixed << std::setprecision(6) << p4_y<< ","<<std::endl;
              }
              row_count = row_count+1;

            }

          }
      }


  }



  void PolyTrajOptimizer::recode_data_cost(const double jerk_cost,const Eigen::VectorXd feas_qvar_costs,const double time_cost){
    if(!is_record_data_){
      return;
    }
    plan_utils::Trajectory traj = jerkOpt_S3NU_container[0].getTraj(1);
    

    std::cout << "Recording data for iteration cost: " << iter_num_ << std::endl;
    std::string file_name = workspace_ + "src/npct_Planner/data_plot/data_raw/data_iter/npct_iter_cost.csv";
    std::cout << "write dir: " << file_name << std::endl;
    if(iter_num_ == 0){std::remove(file_name.c_str());}
    std::ofstream file_out;
    file_out.open(file_name, std::ios::out | std::ios::app);
    if(file_out.is_open()) {
      if(iter_num_ == 0){
        file_out << "iter_num" << ","
          << "jerk_cost" << ","
          << "static_obs_cost" << ","
          << "dynamic_obs_cost" << ","
          << "VAC_cost" << ","
          << "ref_cost" << ","
          << "dotphi_cost" << ","
          << "time_cost" << ","
          << "total_cost" << ","
          <<  "t1" << ","
          <<  "t2" << ","
          <<  "t3" << ","
          <<  "t4" << ","
          <<  "t5" << ","
          <<  "t6" << ","
          <<  "t7" << ","
          <<  "t8" << ","
          <<  "t9" << ","
          <<  "t10" << ","
          <<  "t11" << ","
          <<  "t12" << ","
          << std::endl;
      }
      file_out << iter_num_ << ","
            << std::fixed << std::setprecision(6) << jerk_cost << ","
            << std::fixed << std::setprecision(6) << feas_qvar_costs(0) << ","
            << std::fixed << std::setprecision(6) << feas_qvar_costs(1) << ","
            << std::fixed << std::setprecision(6) << feas_qvar_costs(2)<< ","
            << std::fixed << std::setprecision(6) << feas_qvar_costs(3)<< ","
            << std::fixed << std::setprecision(6) << feas_qvar_costs(4)<< "," 
            << std::fixed << std::setprecision(6) << time_cost<< ","
            << std::fixed << std::setprecision(6) << jerk_cost + feas_qvar_costs.sum() + time_cost<< ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(0) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(1) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(2) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(3) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(4) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(5) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(6) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(7) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(8) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(9) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(10) << ","
            << std::fixed << std::setprecision(6) << traj.getDurations()(11) << ","
            << std::endl;
    }
    file_out.close();



   }

}


