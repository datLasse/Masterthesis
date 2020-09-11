#---------------------------OLD NUMERCIAL METHOD TOBS0 ------------------------------------------







            temp_obs_breakout = 1
            #temp_BB_breakout*coupling_breakout**2*(thompson_coefficent(temp_BB_breakout,density_breakout,time_breakout,time_breakout)**(-2))
            temp_array = np.array([temp_obs_breakout])
        #    print(temp_obs_breakout)
            itcount_breakout = 0
            terminator_breakout = 1

            while terminator_breakout == 1:

                if thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout) == 1:

                    terminator_breakout = 0
                    thompson_coefficent_breakout = 1
                    print("breakout thompson coefficent = 1 | " + str(itcount_breakout)+ " | "+str(temp_obs_breakout))


                else:

                    temp_iter = temp_BB_breakout*coupling_breakout**2*(thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)**(-2))

                    if np.abs(temp_iter-temp_obs_breakout) < 0.0001:

                        temp_obs_breakout = temp_iter
                        thompson_coefficent_breakout = thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)
                        terminator_breakout = 0
                    #    print("breakout solved by iteration | " + str(itcount_breakout)+ " | " +str(temp_obs_breakout))


                    else:
                    #print(thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout))
                        temp_obs_breakout = temp_iter
                        temp_array = np.append(temp_array,temp_obs_breakout)

                        itcount_breakout += 1
        #    print('here it comes')
        #    print(temp_obs_breakout - temp_BB_breakout*coupling_breakout**2*(thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)**(-2)))
        #    print('there it was')
        #    temp_obs_breakout1 = sp.fsolve(temp_break,coupling_breakout**2*temp_BB_breakout,args=(density_breakout,time_breakout,coupling_breakout,temp_BB_breakout))
            temp_obs_breakout2 = sp.root(temp_break,coupling_breakout**2*temp_BB_breakout,args=(density_breakout,time_breakout,coupling_breakout,temp_BB_breakout),method = 'anderson')
            #print('this is the fsolve one'+str(temp_obs_breakout1))
        #    print('this is the root one'+str(temp_obs_breakout2.x))
        #    print('this is mine'+str(temp_obs_breakout))
            x_shu = np.linspace(0,len(temp_array),len(temp_array))
        #    plt.scatter(x_shu, temp_array)
        #    plt.yscale('log')

            plt.clf()


#----------------------------------------------------------------------------------------------------------------------------------------------------------


                terminator = 1
                itcount = 0
                #temp = sp.fsolve(temp_obser, temp_obs,args=(density_breakout,time_breakout,time,coupling_shell,temp_BB_breakout,temp_obs_breakout,thompson_coefficent_breakout))

                #return temp




                while terminator == 1:

                    if thompson_coefficent(temp_obs,density_breakout,time_breakout,time) == thompson_coefficent(temp_obs,density_breakout,time_breakout,time):

                        terminator = 0
                        temp_obs = pre_temp
                        #    print("observable thompson coefficent = 1 | " + str(itcount)+ " | "+str(temp_obs))
                        return temp_obs

                    else:

                        temp_obs_new = temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power(thompson_coefficent(temp_obs,density_breakout,time_breakout,time)/thompson_coefficent_breakout,-2)



                        if np.abs(temp_obs_new-temp_obs) < 0.001:

                            terminator = 0
                            temp_obs = temp_obs_new
                            print("observable solved by iteration | " + str(itcount)+ " | " +str(temp_obs))
                            #x_shu = np.linspace(0,len(temp_array),len(temp_array))
                            #plt.scatter(x_shu, temp_array)
                            #plt.yscale('log')
                            #plt.show()
                            return temp_obs

                        else:
                            #if itcount > 500:
                            #    terminater = 0
                            #    return temp_obs
                            temp_obs = temp_obs_new
                            #print('solving observable' +str(itcount))
                            #temp_array = np.append(temp_array,temp_obs)
                            itcount += 1
