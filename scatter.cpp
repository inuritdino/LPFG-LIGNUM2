#ifndef FSTREAM_H
#include <fstream>
#define FSTREAM_H 1
#endif


int find_branch_by_order(std::vector<Branch> & br, int order, \
			 std::vector<int> & br_ind)
{
  br_ind.resize(0);
  for(int i = 0; i < br.size(); i++){
    if(br[i].order == order)
      br_ind.push_back(i);
  }
  // for(int i = 0; i < br_ind.size(); i++)
  //   std::cout << br_ind[i] << " ";
  // std::cout << std::endl;
  
  return 0;
}

int parent_branch(std::vector<CylData> & cyls, std::vector<Branch> & br, int cc)
{//Find the parent branch of a branch cotaining cyl CC as the base cyl.
  if(br.size() == 0)
    return -1;
  
  //std::cout << "Got here, branche " << br_ind << std::endl;
  //Parent cyl of the 1st cyl of the branch

  int par_cyl = cyls[cc].parent;
  if(par_cyl == -1)
    return -1;

  for(int i = 0; i < br.size(); i++){
    for(int j = 0; j < br[i].cyl_ind.size(); j++){
      if(par_cyl == br[i].cyl_ind[j])
	return i;
    }
  }
  
  return -2;//should not be possible
}

int scatter_output(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{
  std::ofstream scat;
  scat.open("scatter.dat",std::ofstream::out);
  std::vector<int> curr_brs;
  int order = 0;
  float len,rad,ang,gamma,zeta,tot_len;
  int n,m,k;
  V3f rax,newZ,newY,newX,newV;
  V3f newVprojXY,newVprojXZ;

  // Find trunk (branch of 0-order)
  find_branch_by_order(br,order,curr_brs);
  //std::cout << curr_brs.size() << std::endl;

  //*****************************************************************
  //                          HEADER
  //*****************************************************************
  //Write out the header information on Branch and Segment based data-sets
  scat << "# Branch: bra az ltot rini lapar" << std::endl;
  scat << "# Segment: rad len gamma zeta" << std::endl;

  while(curr_brs.size() > 0){
    scat << "# order " << order << std::endl;
    //**********************************************
    //                    BRANCH
    //**********************************************
    if(order > 0){//Define for non-trunk branches
      //scat << "# bra" << std::endl;
      for(int i = 0; i < curr_brs.size(); i++){
	n = br[curr_brs[i]].cyl_ind[0];//1st cyl of the br
	//******* BRA: branching angle of a branch
	ang = (180/M_PI)*acos(cyls[n].axis*cyls[cyls[n].parent].axis);
	if(isnan(ang)){
	  std::cerr << "Bra: instance of NaN: " << std::endl;
	  std::cerr << cyls[n].axis << "; " << cyls[cyls[n].parent].axis <<\
	    ": " << cyls[n].axis*cyls[cyls[n].parent].axis << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << ang << " ";
	}
	//******* AZ: azimuth angle (around parent segment/cylinder)
	//Rotation axis with respect to X, since this direction does not change when
	//transfer the structure to Matlab
	newX = V3f(1.0,0.0,0.0);//global X-axis
	rax = cyls[cyls[n].parent].axis % newX;//rot. axis
	ang = acos(cyls[cyls[n].parent].axis * newX);//ang to rotate
	//global Y-axis coincides with LPFG's (-Z)-axis
	//global Z-axis coincides with LPFG's Y-axis
	//(see also lpfg_to_normal_orientation(V3f&) in methods.cpp)
	newY = V3f(0.0,0.0,-1.0);//global Y-axis = (-Z)-lpfg
	newZ = V3f(0.0,1.0,0.0);//global Z-axis = Y-lpfg
	//The code below is exactly the same (accounting on the inconsistency in the
	//coordinate systems of LPFG and Matlab) as in GEN_SCATTER() function
	//extracting the QSM properties (see Bayes Forest Matlab interface).
	newY = rotate_v3f(newY,rax,-ang);
	newZ = rotate_v3f(newZ,rax,-ang);
	//Approximation, since the real angle would be from projections of
	//cyls[n].axis onto newY-newZ-plane
	ang = (180/M_PI)*atan2(cyls[n].axis*newZ,cyls[n].axis*newY);
	if(isnan(ang)){
	  std::cerr << "Az: instance of NaN: " << std::endl;
	  std::cerr << cyls[n].axis << "; " << cyls[cyls[n].parent].axis <<\
	    "; " << rax << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << ang << " ";
	}
	//******* LTOT: total length of a branch
	tot_len = 0.0;
	for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){
	  tot_len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	}
	if(isnan(tot_len)){
	  std::cerr << "ltot: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << tot_len << " ";
	}
	//******* RINI: Base/inital radius of a branch
	rad = cyls[n].radius;
	if(isnan(rad)){
	  std::cerr << "rini: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << rad << " ";
	}
	//******* LAPAR: Distance from the parent's base (length along the parent)
	m = parent_branch(cyls,br,n);//parent branch of the branch i
	k = cyls[n].parent;//parent cyl of n
	len = 0.0;
	for(int j = 0; j < br[m].cyl_ind.size(); j++){
	  len += cyls[br[m].cyl_ind[j]].length;
	  if(br[m].cyl_ind[j] == k)
	    break;
	}
	if(isnan(len)){
	  std::cerr << "lapar: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << len << " ";
	}

	//******** Trailing newline (NOTE: put new data above this)
	scat << std::endl;
      }
      
      scat << std::endl << std::endl;
    }
    
    //**********************************************
    //                    SEGMENT
    //**********************************************

    for(int i = 0; i < curr_brs.size(); i++){//Over branches
      len = 0.0;
      for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){//Over segments
	//******* RAD: Radius along a branch
	rad = cyls[br[curr_brs[i]].cyl_ind[j]].radius;
	if(isnan(rad)){
	  std::cerr << "RAD: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << rad << " ";
	}
	//******* LEN: length along a branch
	if(isnan(len)){
	  std::cerr << "len: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << len << " ";
	}
	len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	//******* GAMMA and ZETA: horizontal and vertical variations along a branch
	// ** TWO-ROTATION METHOD **
	// 1. Rotate around (Z) axis for X-axis to coincide with the j-1'th cyl(parent)
	// projection onto XY-plane. Rotate the Y-axis correspondingly and
	// do not change Z-axis. Calculate gamma in the newX-newY plane for j'th cyl(child).
	// 2. Rotate around (Y) axis for newX to coincide with the j-1'th cyl. Rotate
	// the Z-axis correspondingly and do not change Y-axis. Calculate zeta in the
	// newX-newZ plane for j'th cyl(child).
	// NOTE: X = X(lpfg), Y = -Z(lpfg), and Z = Y(lpfg)

	if(j == 0){
	  // for the 1st segment make gamma and zeta zeros
	  gamma = 0.0;
	  zeta = 0.0;
	}
	if(j > 0){
	  //For the segments other than the 1st take adjacent pairs
	  //1. Rotation around Z and XY projection angle gamma
	  modify_coord_xy(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,newX,newY);
	  gamma = (180/M_PI) * \
	    projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newY);
	
	  //2. Rotation around Y and XZ projection angle zeta
	  modify_coord_xz(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,\
			  newX,newZ,(const V3f)newY);
	  zeta = (180/M_PI) * \
	    projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newZ);
	}

	if(isnan(gamma)){
	  std::cerr << "gamma: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << gamma << " ";
	}
	if(isnan(zeta)){
	  std::cerr << "zeta: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << zeta << " ";
	}


	//******** Trailing newline (NOTE: put new data above this)
	scat << std::endl;
      }
    }
    scat << std::endl << std::endl;
    
    //+++++++++++ FIND BRANCHES OF THE NEXT ORDER ++++++++++++++
    order++;
    find_branch_by_order(br,order,curr_brs);
  }

  scat.close();
  return 0;
}

int extract_branches(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{
  if(nCyl == 0)
    return 1;

  br.resize(0);
  
  int cc = 0;
  std::vector<int> base_cyls(1,cc);

  Branch br_tmp;
  int n_br = 0;
  
  while(base_cyls.size() > 0){
    cc = base_cyls[0];//Rename the first cyl in BASE_CYLS
    br_tmp.cyl_ind.resize(0);
    br_tmp.cyl_ind.push_back(cc);

    while(cyls[cc].extension > -1){
      for(int i=0;i<cyls[cc].children.size();i++)
	base_cyls.push_back(cyls[cc].children[i]);
      cc = cyls[cc].extension;
      br_tmp.cyl_ind.push_back(cc);
    }
    for(int i=0;i<cyls[cc].children.size();i++)
      base_cyls.push_back(cyls[cc].children[i]);

    base_cyls.erase(base_cyls.begin());

    //Copy to the branch array
    br_tmp.order = cyls[br_tmp.cyl_ind[0]].order;
    br_tmp.parent_br = parent_branch(cyls,br,br_tmp.cyl_ind[0]);
    br.push_back(br_tmp);
    n_br++;
  }

  // for(int i = 0; i < br.size(); i++){
  //   std::cout << "Branch " << i << ":" << std::endl;
  //   std::cout << br[i];
  // }

  return 0;
}
