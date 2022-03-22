function cdf = snow_cdf

cdf_left = [
1.47924078e-31
2.56861807e-30
3.97661847e-29
5.51424866e-28
6.87924067e-27
7.75394565e-26
7.92875882e-25
7.38394964e-24
6.28645449e-23
4.91047138e-22
3.53138261e-21
2.34592381e-20
1.44416193e-19
8.26381067e-19
4.40842497e-18
2.19861921e-17
1.02791130e-16
4.51677693e-16
1.87004463e-15
7.31248580e-15
2.70685782e-14
9.50628202e-14
3.17409400e-13
1.00965908e-12
3.06564732e-12
8.90173424e-12
2.47634940e-11
6.61123493e-11
1.69669992e-10
4.19245871e-10
9.98931368e-10
2.29848036e-09
5.11439043e-09
1.10199603e-08
2.30229438e-08
4.66955370e-08
9.20536345e-08
1.76585117e-07
3.29982895e-07
6.01326509e-07
1.06967017e-06
1.85922268e-06
3.16052586e-06
5.25922645e-06
8.57415904e-06
1.37064794e-05
2.15014735e-05
3.31243868e-05
5.01511559e-05
7.46742667e-05
1.09423133e-04
1.57897413e-04
2.24510608e-04
3.14740191e-04
4.35279459e-04
5.94185417e-04
8.01016242e-04
1.06695153e-03
1.40488841e-03
1.82950683e-03
2.35729812e-03
3.00655182e-03
3.79729698e-03
4.75119597e-03
5.89139019e-03
7.24229910e-03
8.82937548e-03
1.06788217e-02
1.28172732e-02
1.52714557e-02
1.80678258e-02
2.12322014e-02
2.47893934e-02
2.87628449e-02
3.31742877e-02
3.80434238e-02
4.33876368e-02
4.92217416e-02
5.55577734e-02
6.24048222e-02
6.97689130e-02
7.76529317e-02
8.60565982e-02
9.49764834e-02
1.04406068e-01
1.14335842e-01
1.24753435e-01
1.35643781e-01
1.46989309e-01
1.58770152e-01
1.70964379e-01
1.83548228e-01
1.96496362e-01
2.09782111e-01
2.23377732e-01
2.37254647e-01
2.51383691e-01
2.65735337e-01
2.80279915e-01
2.94987822e-01
3.09829708e-01
3.24776652e-01
3.39800320e-01
3.54873109e-01
3.69968268e-01
3.85060006e-01
4.00123587e-01
4.15135402e-01
4.30073032e-01
4.44915294e-01
4.59642275e-01
4.74235353e-01
4.88677208e-01
    ];

cdf_right = [
4.97048180e-01
4.82955537e-01
4.69058316e-01
4.55368718e-01
4.41897724e-01
4.28655129e-01
4.15649586e-01
4.02888647e-01
3.90378807e-01
3.78125559e-01
3.66133440e-01
3.54406083e-01
3.42946267e-01
3.31755976e-01
3.20836441e-01
3.10188199e-01
2.99811137e-01
2.89704546e-01
2.79867165e-01
2.70297228e-01
2.60992507e-01
2.51950356e-01
2.43167752e-01
2.34641330e-01
2.26367423e-01
2.18342097e-01
2.10561180e-01
2.03020294e-01
1.95714888e-01
1.88640258e-01
1.81791576e-01
1.75163913e-01
1.68752261e-01
1.62551550e-01
1.56556669e-01
1.50762484e-01
1.45163848e-01
1.39755621e-01
1.34532680e-01
1.29489930e-01
1.24622315e-01
1.19924829e-01
1.15392521e-01
1.11020504e-01
1.06803964e-01
1.02738160e-01
9.88184365e-02
9.50402195e-02
9.13990264e-02
8.78904661e-02
8.45102422e-02
8.12541549e-02
7.81181024e-02
7.50980822e-02
7.21901916e-02
6.93906281e-02
6.66956895e-02
6.41017739e-02
6.16053786e-02
5.92031000e-02
5.68916324e-02
5.46677670e-02
5.25283907e-02
5.04704845e-02
4.84911223e-02
4.65874689e-02
4.47567789e-02
4.29963942e-02
4.13037429e-02
3.96763369e-02
3.81117701e-02
3.66077166e-02
3.51619288e-02
3.37722349e-02
3.24365376e-02
3.11528116e-02
2.99191019e-02
2.87335218e-02
2.75942506e-02
2.64995323e-02
2.54476731e-02
2.44370399e-02
2.34660582e-02
2.25332104e-02
2.16370339e-02
2.07761196e-02
1.99491098e-02
1.91546970e-02
1.83916215e-02
1.76586706e-02
1.69546768e-02
1.62785157e-02
1.56291054e-02
1.50054045e-02
1.44064106e-02
1.38311595e-02
1.32787233e-02
1.27482092e-02
1.22387586e-02
1.17495456e-02
1.12797758e-02
1.08286853e-02
1.03955393e-02
9.97963163e-03
9.58028302e-03
9.19684056e-03
8.82867656e-03
8.47518763e-03
8.13579378e-03
7.80993757e-03
7.49708323e-03
7.19671585e-03
6.90834064e-03
6.63148211e-03
6.36568336e-03
6.11050538e-03
5.86552638e-03
5.63034110e-03
5.40456020e-03
5.18780964e-03
4.97973011e-03
4.77997645e-03
4.58821709e-03
4.40413355e-03
4.22741992e-03
4.05778237e-03
3.89493871e-03
3.73861788e-03
3.58855958e-03
3.44451381e-03
3.30624048e-03
3.17350902e-03
3.04609802e-03
2.92379484e-03
2.80639532e-03
2.69370339e-03
2.58553078e-03
2.48169671e-03
2.38202760e-03
2.28635678e-03
2.19452422e-03
2.10637627e-03
2.02176538e-03
1.94054991e-03
1.86259388e-03
1.78776669e-03
1.71594301e-03
1.64700246e-03
1.58082951e-03
1.51731322e-03
1.45634710e-03
1.39782890e-03
1.34166047e-03
1.28774758e-03
1.23599977e-03
1.18633019e-03
1.13865549e-03
1.09289563e-03
1.04897378e-03
1.00681619e-03
9.66352077e-04
9.27513464e-04
8.90235116e-04
8.54454408e-04
8.20111224e-04
7.87147855e-04
7.55508907e-04
7.25141207e-04
6.95993712e-04
6.68017427e-04
6.41165321e-04
6.15392249e-04
5.90654877e-04
5.66911609e-04
5.44122518e-04
5.22249279e-04
5.01255104e-04
4.81104682e-04
4.61764118e-04
4.43200879e-04
4.25383736e-04
4.08282715e-04
3.91869044e-04
3.76115107e-04
3.60994396e-04
3.46481467e-04
3.32551900e-04
3.19182253e-04
3.06350026e-04
2.94033622e-04
2.82212314e-04
2.70866204e-04
2.59976196e-04
2.49523959e-04
2.39491900e-04
2.29863132e-04
2.20621445e-04
2.11751282e-04
2.03237711e-04
1.95066400e-04
1.87223590e-04
1.79696080e-04
1.72471194e-04
1.65536769e-04
1.58881129e-04
1.52493068e-04
1.46361830e-04
1.40477092e-04
1.34828944e-04
1.29407876e-04
1.24204759e-04
1.19210832e-04
1.14417685e-04
1.09817247e-04
1.05401772e-04
1.01163822e-04
9.70962619e-05
9.31922411e-05
8.94451852e-05
8.58487837e-05
8.23969800e-05
7.90839607e-05
7.59041463e-05
7.28521815e-05
6.99229262e-05
6.71114471e-05
6.44130090e-05
6.18230672e-05
5.93372597e-05
5.69513996e-05
5.46614687e-05
5.24636102e-05
5.03541220e-05
4.83294515e-05
4.63861882e-05
4.45210593e-05
4.27309232e-05
4.10127647e-05
3.93636900e-05
3.77809213e-05
3.62617928e-05
3.48037457e-05
3.34043241e-05
3.20611708e-05
3.07720235e-05
2.95347108e-05
2.83471485e-05
2.72073363e-05
2.61133542e-05
2.50633596e-05
2.40555839e-05
2.30883294e-05
2.21599670e-05
2.12689329e-05
2.04137262e-05
1.95929063e-05
1.88050906e-05
1.80489521e-05
1.73232171e-05
1.66266631e-05
1.59581168e-05
1.53164521e-05
1.47005881e-05
1.41094874e-05
1.35421543e-05
1.29976331e-05
1.24750066e-05
1.19733945e-05
1.14919517e-05
1.10298674e-05
1.05863631e-05
1.01606918e-05
9.75213632e-06
9.36000858e-06
8.98364799e-06
8.62242058e-06
8.27571786e-06
7.94295581e-06
7.62357388e-06
7.31703408e-06
7.02282004e-06
6.74043615e-06
6.46940673e-06
6.20927523e-06
5.95960345e-06
5.71997081e-06
5.48997366e-06
5.26922455e-06
5.05735163e-06
4.85399799e-06
4.65882109e-06
4.47149214e-06
4.29169558e-06
4.11912853e-06
3.95350032e-06
3.79453192e-06
3.64195555e-06
3.49551420e-06
3.35496118e-06
3.22005972e-06
3.09058258e-06
2.96631165e-06
2.84703759e-06
2.73255948e-06
2.62268448e-06
2.51722750e-06
2.41601089e-06
2.31886415e-06
2.22562363e-06
2.13613228e-06
2.05023932e-06
1.96780008e-06
1.88867569e-06
1.81273285e-06
1.73984363e-06
1.66988526e-06
1.60273988e-06
1.53829439e-06
1.47644022e-06
1.41707318e-06
1.36009326e-06
1.30540448e-06
1.25291471e-06
1.20253553e-06
1.15418207e-06
1.10777289e-06
1.06322980e-06
1.02047777e-06
9.79444776e-07
9.40061704e-07
9.02262208e-07
8.65982614e-07
8.31161807e-07
7.97741129e-07
7.65664283e-07
7.34877234e-07
7.05328118e-07
6.76967161e-07
6.49746585e-07
6.23620537e-07
5.98545007e-07
5.74477754e-07
5.51378235e-07
5.29207538e-07
5.07928316e-07
4.87504722e-07
4.67902353e-07
4.49088186e-07
4.31030530e-07
4.13698965e-07
3.97064295e-07
3.81098498e-07
3.65774680e-07
3.51067026e-07
3.36950761e-07
3.23402105e-07
3.10398235e-07
2.97917245e-07
2.85938111e-07
2.74440653e-07
2.63405503e-07
2.52814073e-07
2.42648519e-07
2.32891718e-07
2.23527235e-07
2.14539293e-07
2.05912753e-07
1.97633083e-07
1.89686335e-07
1.82059122e-07
1.74738597e-07
1.67712427e-07
1.60968776e-07
1.54496285e-07
1.48284051e-07
1.42321607e-07
1.36598912e-07
1.31106323e-07
1.25834590e-07
1.20774831e-07
1.15918522e-07
1.11257484e-07
1.06783865e-07
1.02490127e-07
9.83690397e-08
9.44136593e-08
9.06173232e-08
8.69736363e-08
8.34764606e-08
8.01199049e-08
7.68983150e-08
7.38062641e-08
7.08385431e-08
6.79901532e-08
6.52562959e-08
6.26323658e-08
6.01139430e-08
5.76967850e-08
5.53768199e-08
5.31501397e-08
5.10129935e-08
4.89617811e-08
4.69930471e-08
4.51034751e-08
4.32898821e-08
4.15492130e-08
3.98785354e-08
3.82750351e-08
3.67360110e-08
3.52588704e-08
3.38411250e-08
3.24803867e-08
3.11743631e-08
2.99208541e-08
2.87177484e-08
2.75630190e-08
2.64547207e-08
2.53909869e-08
2.43700252e-08
2.33901161e-08
2.24496086e-08
2.15469187e-08
2.06805255e-08
1.98489697e-08
1.90508505e-08
1.82848232e-08
1.75495977e-08
1.68439351e-08
1.61666471e-08
1.55165926e-08
1.48926765e-08
1.42938478e-08
1.37190977e-08
1.31674582e-08
1.26379999e-08
1.21298309e-08
1.16420952e-08
1.11739711e-08
1.07246700e-08
1.02934353e-08
9.87954030e-09
9.48228784e-09
9.10100872e-09
8.73506079e-09
8.38382741e-09
8.04671696e-09
7.72316167e-09
7.41261630e-09
7.11455794e-09
6.82848433e-09
6.55391363e-09
6.29038333e-09
6.03744954e-09
5.79468606e-09
5.56168400e-09
5.33805089e-09
5.12340992e-09
4.91739960e-09
4.71967287e-09
4.52989668e-09
4.34775138e-09
4.17293000e-09
4.00513811e-09
3.84409304e-09
3.68952358e-09
3.54116925e-09
3.39878026e-09
3.26211658e-09
3.13094817e-09
3.00505398e-09
2.88422197e-09
2.76824852e-09
2.65693834e-09
2.55010391e-09
2.44756515e-09
2.34914954e-09
2.25469121e-09
2.16403095e-09
2.07701611e-09
1.99350014e-09
1.91334226e-09
1.83640747e-09
1.76256632e-09
1.69169412e-09
1.62367175e-09
1.55838464e-09
1.49572255e-09
1.43558010e-09
1.37785605e-09
1.32245292e-09
1.26927768e-09
1.21824051e-09
1.16925558e-09
1.12224019e-09
1.07711540e-09
1.03380504e-09
9.92236071e-10
9.52338652e-10
9.14045506e-10
8.77292128e-10
8.42016568e-10
8.08159428e-10
7.75663644e-10
7.44474482e-10
7.14539428e-10
6.85808077e-10
6.58232024e-10
6.31764752e-10
6.06361739e-10
5.81980131e-10
5.58578961e-10
5.36118705e-10
5.14561616e-10
4.93871277e-10
4.74012940e-10
4.54953075e-10
4.36659597e-10
4.19101642e-10
4.02249789e-10
3.86075394e-10
3.70551478e-10
3.55651730e-10
3.41351170e-10
3.27625593e-10
3.14451909e-10
3.01807912e-10
2.89672286e-10
2.78024714e-10
2.66845435e-10
2.56115684e-10
2.45817366e-10
2.35933162e-10
2.26446417e-10
2.17341034e-10
2.08601802e-10
2.00214068e-10
1.92163507e-10
1.84436688e-10
1.77020620e-10
1.69902648e-10
1.63070890e-10
1.56513913e-10
1.50220503e-10
1.44180223e-10
1.38382750e-10
1.32818534e-10
1.27477917e-10
1.22352017e-10
1.17432286e-10
1.12710397e-10
1.08178355e-10
1.03828612e-10
9.96536187e-11
9.56466017e-11
9.18006782e-11
8.81094087e-11
8.45665760e-11
8.11661849e-11
7.79025733e-11
7.47700790e-11
7.17635951e-11
6.88780144e-11
6.61084520e-11
6.34503561e-11
6.08989525e-11
5.84502446e-11
5.61000135e-11
5.38442624e-11
5.16792165e-11
4.96012120e-11
4.76066964e-11
4.56924498e-11
4.38552528e-11
4.20917745e-11
4.03993505e-11
3.87748722e-11
3.72157860e-11
3.57193164e-11
3.42830209e-11
3.29045680e-11
3.15815152e-11
3.03116421e-11
2.90927282e-11
2.79229972e-11
2.68002287e-11
2.57225352e-11
2.46882514e-11
2.36956010e-11
2.27428076e-11
2.18283169e-11
2.09505746e-11
2.01081374e-11
1.92996730e-11
1.85236271e-11
1.77787784e-11
1.70639058e-11
1.63777880e-11
1.57192037e-11
1.50871537e-11
1.44805279e-11
1.38982159e-11
1.33394407e-11
1.28030919e-11
1.22882815e-11
1.17941212e-11
1.13199450e-11
1.08647535e-11
1.04278808e-11
1.00085495e-11
9.60609370e-12
9.21984711e-12
8.84914364e-12
8.49331716e-12
8.15181256e-12
7.82407472e-12
7.50943752e-12
7.20745685e-12
6.91768864e-12
6.63946675e-12
6.37256914e-12
6.11632966e-12
5.87041527e-12
5.63427083e-12
5.40778533e-12
5.19029264e-12
4.98157071e-12
4.78128648e-12
4.58910687e-12
4.40447678e-12
4.22739621e-12
4.05742107e-12
3.89432930e-12
3.73767683e-12
3.58746366e-12
3.44313467e-12
3.30468986e-12
3.17179616e-12
3.04434256e-12
2.92188496e-12
2.80442336e-12
2.69162470e-12
2.58337796e-12
2.47957210e-12
2.37987408e-12
2.28417285e-12
2.19224638e-12
2.10409468e-12
2.01949568e-12
1.93833838e-12
1.86040072e-12
1.78557169e-12
1.71374026e-12
1.64490643e-12
1.57873714e-12
1.51523238e-12
1.45428114e-12
1.39588341e-12
1.33970612e-12
1.28586031e-12
1.23412391e-12
1.18449694e-12
1.13686838e-12
1.09112719e-12
1.04727338e-12
1.00519593e-12
9.64783808e-13
9.26037025e-13
8.88733531e-13
8.52984350e-13
8.18678458e-13
7.85815857e-13
7.54174501e-13
7.23865412e-13
6.94777569e-13
6.66799949e-13
6.40043574e-13
6.14286400e-13
5.89528426e-13
5.65880676e-13
5.43121104e-13
5.21249710e-13
5.00266495e-13
4.80171458e-13
4.60853578e-13
4.42312853e-13
4.24549285e-13
4.07451850e-13
3.91131572e-13
3.75366405e-13
3.60267371e-13
3.45834472e-13
3.31956684e-13
3.18522986e-13
3.05755421e-13
2.93431945e-13
2.81663581e-13
2.70339306e-13
2.59459121e-13
2.49023024e-13
2.39031017e-13
2.29372077e-13
2.20157226e-13
2.11275442e-13
2.02837747e-13
1.94622096e-13
1.86850535e-13
1.79301018e-13
1.72084569e-13
1.65201186e-13
1.58539848e-13
1.52211577e-13
1.46105350e-13
1.40221168e-13
1.34559031e-13
1.29118938e-13
1.23900890e-13
1.19015908e-13
1.14130927e-13
1.09579013e-13
1.05138120e-13
1.00919273e-13
9.69224700e-14
9.30366895e-14
8.92619312e-14
8.57092175e-14
8.22675261e-14
7.89368571e-14
7.57172103e-14
7.27196081e-14
6.98330282e-14
6.69464484e-14
6.42819131e-14
6.17284002e-14
5.91748872e-14
5.68434189e-14
5.45119505e-14
5.24025268e-14
5.02931030e-14
4.81836793e-14
4.62963001e-14
4.44089210e-14
4.26325641e-14
4.09672296e-14
3.93018951e-14
3.77475828e-14
3.61932706e-14
3.47499807e-14
3.33066907e-14
3.19744231e-14
3.07531778e-14
2.94209102e-14
2.83106871e-14
2.72004641e-14
2.60902411e-14
2.49800181e-14
2.39808173e-14
2.30926389e-14
2.20934382e-14
2.12052598e-14
2.03170814e-14
1.95399252e-14
1.87627691e-14
1.79856130e-14
1.73194792e-14
1.65423231e-14
1.58761893e-14
1.53210777e-14
1.46549439e-14
1.40998324e-14
1.35447209e-14
1.29896094e-14
1.24344979e-14
1.19904087e-14
1.14352972e-14
1.09912079e-14
1.05471187e-14
1.01030295e-14
9.76996262e-15
9.32587341e-15
8.99280650e-15
8.65973959e-15
8.21565038e-15
7.88258347e-15
7.66053887e-15
7.32747196e-15
6.99440506e-15
6.77236045e-15
6.43929354e-15
6.21724894e-15
5.99520433e-15
5.66213743e-15
5.44009282e-15
5.21804822e-15
4.99600361e-15
4.88498131e-15
4.66293670e-15
4.44089210e-15
4.32986980e-15
4.10782519e-15
3.99680289e-15
3.77475828e-15
3.66373598e-15
3.44169138e-15
3.33066907e-15
3.21964677e-15
3.10862447e-15
2.99760217e-15
2.88657986e-15
2.77555756e-15
2.66453526e-15
2.55351296e-15
2.44249065e-15
2.33146835e-15
2.22044605e-15
2.10942375e-15
1.99840144e-15
1.88737914e-15
1.77635684e-15
1.66533454e-15
1.55431223e-15
1.44328993e-15
1.33226763e-15
1.22124533e-15
1.11022302e-15
9.99200722e-16
8.88178420e-16
7.77156117e-16
6.66133815e-16
5.55111512e-16
4.44089210e-16
3.33066907e-16
2.22044605e-16
1.11022302e-16
    ];

cdf = [cdf_left; 1-cdf_right];

end