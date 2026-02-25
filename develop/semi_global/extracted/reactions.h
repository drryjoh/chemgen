    
double call_forward_reaction_0(const double& temperature, const double& log_temperature)  
{ 
	return 
	arrhenius(double(1.0000000000000002e-10), 
			  double(0.0), 
			  double(41840000.0), 
			  temperature, 
			  log_temperature);}
    
double call_forward_reaction_1(const double& temperature, const double& log_temperature)  
{ 
	return 
	arrhenius(double(200000000.00000003), 
			  double(0.0), 
			  double(83680000.0), 
			  temperature, 
			  log_temperature);
}
    
double call_forward_reaction_2(const double& temperature, const double& log_temperature)  
{ 
	return 
	arrhenius(double(300000000000.00006), 
			  double(0.0), 
			  double(125520000.0), 
			  temperature, 
			  log_temperature);
}

#if 0
double dcall_forward_reaction_0_dtemperature(const double& temperature, const double& log_temperature)  { return darrhenius_dtemperature(double(1.0000000000000002e-10), double(0.0), double(41840000.0), temperature, log_temperature);}
    
double dcall_forward_reaction_0_dlog_temperature(const double& temperature, const double& log_temperature)  { return darrhenius_dlog_temperature(double(1.0000000000000002e-10), double(0.0), double(41840000.0), temperature, log_temperature);}
    
double dcall_forward_reaction_1_dtemperature(const double& temperature, const double& log_temperature)  { return darrhenius_dtemperature(double(200000000.00000003), double(0.0), double(83680000.0), temperature, log_temperature);}
    
double dcall_forward_reaction_1_dlog_temperature(const double& temperature, const double& log_temperature)  { return darrhenius_dlog_temperature(double(200000000.00000003), double(0.0), double(83680000.0), temperature, log_temperature);}
    
double dcall_forward_reaction_2_dtemperature(const double& temperature, const double& log_temperature)  { return darrhenius_dtemperature(double(300000000000.00006), double(0.0), double(125520000.0), temperature, log_temperature);}
    
double dcall_forward_reaction_2_dlog_temperature(const double& temperature, const double& log_temperature)  { return darrhenius_dlog_temperature(double(300000000000.00006), double(0.0), double(125520000.0), temperature, log_temperature);}
#endif
