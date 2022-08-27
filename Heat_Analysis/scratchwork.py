# cp1,_,_ =heat_analysis.properties_He(6.8,Constants.SET_PRESSURE_MOCHI)
# print(cp1)
# cp2,_,_ =heat_analysis.properties_He(40,Constants.SET_PRESSURE_MOCHI)
# print(cp2)
# cps,_,_ =heat_analysis.properties_He(20,Constants.SET_PRESSURE_MOCHI)
# print(cps)
# cpl,_,_ =heat_analysis.properties_He(60,Constants.SET_PRESSURE_MOCHI)
# print(cpl)
# cpu,_,_ =heat_analysis.properties_He(282,Constants.SET_PRESSURE_MOCHI)
# print(cpu)
# cpx,_,_ =heat_analysis.properties_He(133,Constants.SET_PRESSURE_MOCHI)
# print(cpx)
# cp_ave=((cp1+cp2)/2)
# cp_ave=(cps+cp2)/2
# cp_ave=(cpl+cps)/2
# cp_ave=(cpu+cpl)/2
# cp_ave=(cpx+cpu)/2
# cp_ave=(cpx+cp1)/2
# print(cp_ave)
# f=0.33
# print(cp_ave*m_dot()*(40-133))
# Tr=f*282+(1-f)*60
# # print(Tr)



# convective cooling = h*A*(T_surface - T_fluid)