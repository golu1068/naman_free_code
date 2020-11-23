import numpy as np
###########################################
def antoine_eqn(T):
    global A, B, C
    A = 8.07131;
    B = 1730.63;
    C = 233.426;
    psat_mmHg = 10**(A-B/(T+C)); #mmHg
    pwv_sat = psat_mmHg/760*101325; #Pa
    return pwv_sat

def heat1_bc(xl,ul,xr,ur,t):
    pl = h*(T_inf - ul);
    ql = rho*cp; ##K DT/Dx
    K = alpha*rho*cp
    pr = -h*(T_inf - ur);
    qr = rho*cp;
    return pl,ql,pr,qr

def heat1_ic(x):
    u0 = 25 + 273  ##degree celsius
    return u0

def heat1_pde(x,t,u,DuDx):
    c = 1
    f = K/(rho*cp)*DuDx
    s = 0
    return c,f,s

def hheqn(u):
    A = 1.2733890;
    B = 0.121262;
    C = 0.001009;
    phi = ((B*u-1)+np.sqrt((1-B*u)**2+4*A*C*u**2))/(2*C*u)
    return phi

def mass1_bc(xl,ul,xr,ur,t):
    pl = -hm*(Pwv_inf - RH*Pwv_sat)
    ql = rho
    pr = hm*(Pwv_inf - RH*Pwv_sat)
    qr = rho
    return pl,ql,pr,qr

def mass1_ic(x):
    u0 = 0.25
    return u0

def mass1_pde(x,t,u,DuDx):
    c = 1;
    f = Deff*DuDx;
    s = 0;
    return c,f,s

def MC_conversion(m_water,m_solid):
    n = length(m_water)
    MC_wb=[];MC_db=[];
    for i in range(0, n):
        MC_wb[i] = m_water[i]/(m_solid[i]+m_water[i])
        MC_db[i] = m_water[i]/m_solid[i]
    return MC_wb, MC_db
######################################################
moist_sol = [[0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25], [0.249037577394107, 0.249502961729266, 0.249764874876708, 0.249895419203553, 0.249955381535084, 0.249981480633068, 0.249992443346379, 0.249996931472359, 0.249998709660432, 0.249999325309582, 0.249999325309582, 0.249998709660432, 0.249996931472359, 0.249992443346379, 0.249981480633068, 0.249955381535084, 0.249895419203553, 0.249764874876708, 0.249502961729266, 0.249037577394108], [0.248602952527691, 0.249106919837289, 0.249460767881262, 0.249690878577939, 0.249830606441323, 0.249910627527552, 0.249954221129253, 0.249976878057693, 0.249987910146879, 0.249992347063935, 0.249992347063935, 0.249987910146879, 0.249976878057693, 0.249954221129253, 0.249910627527552, 0.249830606441323, 0.249690878577939, 0.249460767881262, 0.249106919837289, 0.248602952527691], [0.248308035537278, 0.248823854675883, 0.249210067439713, 0.249484510034287, 0.249670879769694, 0.24979270896295, 0.249869659649688, 0.249916334221883, 0.249942583905046, 0.249954332962674, 0.249954332962674, 0.249942583905046, 0.249916334221883, 0.249869659649688, 0.24979270896295, 0.249670879769694, 0.249484510034287, 0.249210067439713, 0.248823854675883, 0.248308035537278], [0.248020691749524, 0.248547531099927, 0.24896393147318, 0.249280023996246, 0.249510646779966, 0.249672592305114, 0.24978192866972, 0.249852184824937, 0.249893518472849, 0.24991256235724, 0.24991256235724, 0.249893518472849, 0.249852184824937, 0.24978192866972, 0.249672592305115, 0.249510646779966, 0.249280023996246, 0.24896393147318, 0.248547531099927, 0.248020691749524], [0.247785221931974, 0.24831770788079, 0.248750082037152, 0.249089822469135, 0.249348081012035, 0.249537871639471, 0.249672159415152, 0.249762361804868, 0.24981743543309, 0.249843449659747, 0.249843449659747, 0.24981743543309, 0.249762361804868, 0.249672159415152, 0.249537871639471, 0.249348081012035, 0.249089822469135, 0.248750082037152, 0.24831770788079, 0.247785221931974], [0.247557563511438, 0.248095114322378, 0.248541849698801, 0.248902949047007, 0.249186338380517, 0.249401663056548, 0.249559039507227, 0.24966786005922, 0.2497358518749, 0.249768451078597, 0.249768451078597, 0.2497358518749, 0.24966786005922, 0.249559039507227, 0.249401663056548, 0.249186338380517, 0.248902949047006, 0.248541849698801, 0.248095114322378, 0.247557563511438], [0.247355254250333, 0.247895982101549, 0.248351845688146, 0.248726875826638, 0.249027266950599, 0.24926062595802, 0.249435046224668, 0.249558175144617, 0.249636418306218, 0.249674351918777, 0.249674351918777, 0.249636418306218, 0.249558175144617, 0.249435046224668, 0.24926062595802, 0.249027266950599, 0.248726875826638, 0.248351845688146, 0.247895982101549, 0.247355254250333], [0.247164179843823, 0.247707422196146, 0.248170517527737, 0.248556606029941, 0.248870534685548, 0.249118311193242, 0.249306419495971, 0.249441091746506, 0.24952763703337, 0.249569901877512, 0.249569901877512, 0.24952763703337, 0.249441091746506, 0.249306419495971, 0.249118311193242, 0.248870534685548, 0.248556606029941, 0.248170517527737, 0.247707422196146, 0.247164179843823], [0.246993242823555, 0.247537812139835, 0.248004739986686, 0.248396738306806, 0.248717995146449, 0.24897370633574, 0.249169487764862, 0.249310747284915, 0.249402100908232, 0.249446898886068, 0.249446898886068, 0.249402100908232, 0.249310747284915, 0.249169487764862, 0.24897370633574, 0.248717995146449, 0.248396738306806, 0.248004739986686, 0.247537812139835, 0.246993242823555], [0.246822305803287, 0.247368202083523, 0.247838962445634, 0.248236870583672, 0.248565455607349, 0.248829101478239, 0.249032556033753, 0.249180402823323, 0.249276564783094, 0.249323895894624, 0.249323895894624, 0.249276564783094, 0.249180402823323, 0.249032556033753, 0.248829101478239, 0.248565455607349, 0.248236870583672, 0.247838962445634, 0.247368202083523, 0.246822305803287], [0.246651368783019, 0.247198592027212, 0.247673184904583, 0.248077002860537, 0.24841291606825, 0.248684496620737, 0.248895624302644, 0.249050058361732, 0.249151028657956, 0.249200892903181, 0.249200892903181, 0.249151028657956, 0.249050058361732, 0.248895624302644, 0.248684496620737, 0.24841291606825, 0.248077002860537, 0.247673184904583, 0.247198592027212, 0.246651368783019], [0.246480431762751, 0.247028981970901, 0.247507407363532, 0.247917135137402, 0.248260376529151, 0.248539891763235, 0.248758692571535, 0.24891971390014, 0.249025492532818, 0.249077889911737, 0.249077889911737, 0.249025492532818, 0.24891971390014, 0.248758692571535, 0.248539891763235, 0.248260376529151, 0.247917135137402, 0.247507407363532, 0.247028981970901, 0.246480431762751], [0.246316308309803, 0.246865809072326, 0.247346980995184, 0.247760947560942, 0.248109451760383, 0.248394674722418, 0.248619005355129, 0.248784782681926, 0.248894038333108, 0.248948269388995, 0.248948269388995, 0.248894038333108, 0.248784782681926, 0.248619005355129, 0.248394674722418, 0.248109451760383, 0.247760947560942, 0.247346980995184, 0.246865809072326, 0.246316308309803], [0.246160244355546, 0.246710250432469, 0.247192884316904, 0.2476091130841, 0.24796043703943, 0.248248733554222, 0.248476058784392, 0.248644425971124, 0.248755583877582, 0.248810821250802, 0.248810821250802, 0.248755583877582, 0.248644425971124, 0.248476058784392, 0.248248733554222, 0.24796043703943, 0.2476091130841, 0.247192884316904, 0.246710250432469, 0.246160244355546], [0.24600418040129, 0.246554691792611, 0.247038787638624, 0.247457278607259, 0.247811422318477, 0.248102792386027, 0.248333112213655, 0.248504069260321, 0.248617129422056, 0.248673373112608, 0.248673373112608, 0.248617129422056, 0.248504069260321, 0.248333112213655, 0.248102792386027, 0.247811422318477, 0.247457278607259, 0.247038787638624, 0.246554691792611, 0.24600418040129], [0.245848116447034, 0.246399133152754, 0.246884690960345, 0.247305444130417, 0.247662407597524, 0.247956851217832, 0.248190165642918, 0.248363712549519, 0.24847867496653, 0.248535924974414, 0.248535924974414, 0.24847867496653, 0.248363712549519, 0.248190165642918, 0.247956851217832, 0.247662407597524, 0.247305444130417, 0.246884690960345, 0.246399133152754, 0.245848116447034], [0.245692052492777, 0.246243574512896, 0.246730594282065, 0.247153609653576, 0.247513392876571, 0.247810910049636, 0.248047219072181, 0.248223355838716, 0.248340220511004, 0.248398476836221, 0.248398476836221, 0.248340220511004, 0.248223355838716, 0.248047219072181, 0.247810910049636, 0.247513392876571, 0.247153609653576, 0.246730594282065, 0.246243574512896, 0.245692052492777], [0.245539309994701, 0.246091157912873, 0.246579120598455, 0.247003595130786, 0.24736519736092, 0.247664697570914, 0.247902939050262, 0.248080747426341, 0.248198840122672, 0.248257746172167, 0.248257746172167, 0.248198840122672, 0.248080747426341, 0.247902939050262, 0.247664697570914, 0.24736519736092, 0.247003595130786, 0.246579120598455, 0.246091157912873, 0.245539309994701], [0.245389197249617, 0.245941229013488, 0.246429723663028, 0.24685502155117, 0.247217650448617, 0.247518270282892, 0.247757603272358, 0.247936356235984, 0.248055143135796, 0.248114416578086, 0.248114416578086, 0.248055143135796, 0.247936356235984, 0.247757603272358, 0.247518270282892, 0.247217650448617, 0.24685502155117, 0.246429723663028, 0.245941229013488, 0.245389197249617], [0.245239084504532, 0.245791300114102, 0.246280326727601, 0.246706447971553, 0.247070103536314, 0.247371842994871, 0.247612267494454, 0.247791965045626, 0.24791144614892, 0.247971086984004, 0.247971086984004, 0.24791144614892, 0.247791965045626, 0.247612267494454, 0.247371842994871, 0.247070103536314, 0.246706447971553, 0.246280326727601, 0.245791300114102, 0.245239084504532], [0.245088971759448, 0.245641371214716, 0.246130929792173, 0.246557874391937, 0.246922556624011, 0.24722541570685, 0.24746693171655, 0.247647573855268, 0.247767749162045, 0.247827757389923, 0.247827757389923, 0.247767749162045, 0.247647573855268, 0.24746693171655, 0.24722541570685, 0.246922556624011, 0.246557874391937, 0.246130929792173, 0.245641371214716, 0.245088971759448], [0.244938859014364, 0.24549144231533, 0.245981532856746, 0.246409300812321, 0.246775009711708, 0.247078988418829, 0.247321595938647, 0.247503182664911, 0.247624052175169, 0.247684427795841, 0.247684427795841, 0.247624052175169, 0.247503182664911, 0.247321595938647, 0.247078988418829, 0.246775009711708, 0.246409300812321, 0.245981532856746, 0.24549144231533, 0.244938859014364], [0.244790200365971, 0.245342889206804, 0.245833285088894, 0.246261525512691, 0.246627823322257, 0.246932443967087, 0.247175677002716, 0.247357804802419, 0.247479071816556, 0.247539657757759, 0.247539657757759, 0.247479071816556, 0.247357804802419, 0.247175677002716, 0.246932443967087, 0.246627823322257, 0.246261525512691, 0.245833285088894, 0.245342889206804, 0.24479020036597], [0.244642297109772, 0.245195050811188, 0.245685634304876, 0.246114164913455, 0.246480824221692, 0.246785838649683, 0.247029455120617, 0.247211914371241, 0.247333424756066, 0.247394139419951, 0.247394139419951, 0.247333424756066, 0.247211914371241, 0.247029455120617, 0.246785838649683, 0.246480824221691, 0.246114164913455, 0.245685634304876, 0.245195050811188, 0.244642297109772], [0.244494393853573, 0.245047212415572, 0.245537983520858, 0.245966804314219, 0.246333825121126, 0.246639233332278, 0.246883233238518, 0.247066023940064, 0.247187777695576, 0.247248621082143, 0.247248621082143, 0.247187777695576, 0.247066023940064, 0.246883233238518, 0.246639233332278, 0.246333825121126, 0.245966804314219, 0.245537983520858, 0.245047212415572, 0.244494393853573], [0.244346490597374, 0.244899374019955, 0.245390332736839, 0.245819443714983, 0.246186826020561, 0.246492628014873, 0.246737011356419, 0.246920133508887, 0.247042130635086, 0.247103102744335, 0.247103102744335, 0.247042130635086, 0.246920133508887, 0.246737011356419, 0.246492628014873, 0.246186826020561, 0.245819443714983, 0.245390332736839, 0.244899374019955, 0.244346490597374], [0.244198587341175, 0.244751535624339, 0.245242681952821, 0.245672083115747, 0.246039826919995, 0.246346022697468, 0.246590789474319, 0.246774243077709, 0.246896483574597, 0.246957584406527, 0.246957584406527, 0.246896483574597, 0.24677424307771, 0.246590789474319, 0.246346022697468, 0.246039826919995, 0.245672083115747, 0.245242681952821, 0.244751535624339, 0.244198587341175], [0.244051280458885, 0.244604261394206, 0.245095502156792, 0.24552504934405, 0.245892974989564, 0.24619936875894, 0.24644432823348, 0.246627948335334, 0.246750311040635, 0.246811476495153, 0.246811476495153, 0.246750311040635, 0.246627948335334, 0.24644432823348, 0.24619936875894, 0.245892974989564, 0.24552504934405, 0.245095502156792, 0.244604261394206, 0.244051280458885], [0.243904163856815, 0.244457167167812, 0.244948472635109, 0.245378119850583, 0.245746170015523, 0.246052699307261, 0.246297790622375, 0.24648152459264, 0.246603970848084, 0.246665180473289, 0.246665180473289, 0.246603970848084, 0.24648152459264, 0.246297790622375, 0.246052699307261, 0.245746170015523, 0.245378119850583, 0.244948472635109, 0.244457167167812, 0.243904163856815], [0.243757047254745, 0.244310072941417, 0.244801443113425, 0.245231190357117, 0.245599365041481, 0.245906029855582, 0.246151253011269, 0.246335100849946, 0.246457630655533, 0.246518884451424, 0.246518884451424, 0.246457630655533, 0.246335100849946, 0.246151253011269, 0.245906029855582, 0.245599365041481, 0.245231190357117, 0.244801443113425, 0.244310072941417, 0.243757047254745], [0.243609930652675, 0.244162978715023, 0.244654413591742, 0.24508426086365, 0.24545256006744, 0.245759360403903, 0.246004715400164, 0.246188677107252, 0.246311290462982, 0.24637258842956, 0.24637258842956, 0.246311290462982, 0.246188677107252, 0.246004715400164, 0.245759360403903, 0.24545256006744, 0.24508426086365, 0.244654413591742, 0.244162978715023, 0.243609930652675], [0.243462814050606, 0.244015884488629, 0.244507384070058, 0.244937331370183, 0.245305755093398, 0.245612690952225, 0.245858177789059, 0.246042253364558, 0.246164950270431, 0.246226292407695, 0.246226292407695, 0.246164950270431, 0.246042253364558, 0.245858177789059, 0.245612690952225, 0.245305755093398, 0.244937331370184, 0.244507384070058, 0.244015884488629, 0.243462814050606], [0.243315931587417, 0.243869011722415, 0.244360539343189, 0.244790529984531, 0.245159007651677, 0.245466002209565, 0.245711546139108, 0.24589567101523, 0.246018404089506, 0.246079765346094, 0.246079765346094, 0.246018404089506, 0.24589567101523, 0.245711546139108, 0.245466002209565, 0.245159007651677, 0.244790529984531, 0.244360539343189, 0.243869011722415, 0.243315931587417], [0.243169087838072, 0.243722175573676, 0.244213725171336, 0.244643749780948, 0.245012269722674, 0.245319310277226, 0.245564898940245, 0.245749062440987, 0.245871823849302, 0.24593320008308, 0.24593320008308, 0.245871823849302, 0.245749062440987, 0.245564898940245, 0.245319310277226, 0.245012269722675, 0.244643749780948, 0.244213725171336, 0.243722175573676, 0.243169087838072], [0.243022244088727, 0.243575339424937, 0.244066910999483, 0.244496969577364, 0.244865531793672, 0.245172618344886, 0.245418251741382, 0.245602453866744, 0.245725243609097, 0.245786634820066, 0.245786634820066, 0.245725243609097, 0.245602453866744, 0.245418251741382, 0.245172618344886, 0.244865531793672, 0.244496969577364, 0.244066910999483, 0.243575339424937, 0.243022244088727], [0.242875400339382, 0.243428503276198, 0.24392009682763, 0.244350189373781, 0.244718793864669, 0.245025926412547, 0.245271604542519, 0.245455845292501, 0.245578663368893, 0.245640069557052, 0.245640069557052, 0.245578663368893, 0.245455845292501, 0.245271604542519, 0.245025926412547, 0.244718793864669, 0.244350189373781, 0.24392009682763, 0.243428503276198, 0.242875400339382], [0.242728556590037, 0.243281667127459, 0.243773282655777, 0.244203409170197, 0.244572055935667, 0.244879234480208, 0.245124957343656, 0.245309236718258, 0.245432083128688, 0.245493504294039, 0.245493504294039, 0.245432083128688, 0.245309236718258, 0.245124957343656, 0.244879234480208, 0.244572055935667, 0.244203409170197, 0.243773282655777, 0.243281667127459, 0.242728556590037], [0.242581802001248, 0.243134915306147, 0.243626538836585, 0.244056677719563, 0.244425339878433, 0.244732535172152, 0.244978274325291, 0.245162567765452, 0.245285424494268, 0.245346851113938, 0.245346851113938, 0.245285424494268, 0.245162567765452, 0.244978274325291, 0.244732535172152, 0.244425339878433, 0.244056677719563, 0.243626538836585, 0.243134915306147, 0.242581802001248], [0.242435051310254, 0.242988167171341, 0.243479798092972, 0.24390994840024, 0.244278624777358, 0.244585835541655, 0.24483158974102, 0.245015896173102, 0.24513876243272, 0.245200194090402, 0.245200194090402, 0.24513876243272, 0.245015896173102, 0.24483158974102, 0.244585835541655, 0.244278624777358, 0.243909948400241, 0.243479798092972, 0.242988167171341, 0.242435051310254], [0.24228830061926, 0.242841419036535, 0.243333057349358, 0.243763219080918, 0.244131909676283, 0.244439135911159, 0.244684905156749, 0.244869224580752, 0.244992100371173, 0.245053537066867, 0.245053537066867, 0.244992100371173, 0.244869224580752, 0.244684905156749, 0.244439135911159, 0.244131909676283, 0.243763219080918, 0.243333057349358, 0.242841419036535, 0.24228830061926], [0.242141549928266, 0.242694670901729, 0.243186316605744, 0.243616489761596, 0.243985194575208, 0.244292436280662, 0.244538220572478, 0.244722552988401, 0.244845438309625, 0.244906880043331, 0.244906880043331, 0.244845438309625, 0.244722552988401, 0.244538220572478, 0.244292436280662, 0.243985194575208, 0.243616489761596, 0.243186316605744, 0.242694670901729, 0.242141549928266], [0.241994801061733, 0.242547924492483, 0.243039577301728, 0.243469761439884, 0.243838479921684, 0.244145736499237, 0.244391535255248, 0.24457588016055, 0.24469877464393, 0.244760221220785, 0.244760221220785, 0.24469877464393, 0.24457588016055, 0.244391535255248, 0.244145736499237, 0.243838479921684, 0.243469761439884, 0.243039577301728, 0.242547924492483, 0.241994801061733], [0.241848081764965, 0.242401206050097, 0.242892861329855, 0.243323049286828, 0.243691772521787, 0.243999034271665, 0.244244838058643, 0.244429187308434, 0.244552084979155, 0.244613533240934, 0.244613533240934, 0.244552084979155, 0.244429187308435, 0.244244838058643, 0.243999034271665, 0.243691772521787, 0.243323049286828, 0.242892861329855, 0.242401206050097, 0.241848081764965], [0.241701362468198, 0.24225448760771, 0.242746145357983, 0.243176337133772, 0.24354506512189, 0.243852332044093, 0.244098140862037, 0.244282494456319, 0.24440539531438, 0.244466845261082, 0.244466845261082, 0.24440539531438, 0.244282494456319, 0.244098140862037, 0.243852332044093, 0.24354506512189, 0.243176337133772, 0.242746145357983, 0.24225448760771, 0.241701362468198], [0.24155464317143, 0.242107769165324, 0.242599429386111, 0.243029624980717, 0.243398357721994, 0.243705629816521, 0.243951443665432, 0.244135801604203, 0.244258705649605, 0.244320157281231, 0.244320157281231, 0.244258705649605, 0.244135801604203, 0.243951443665432, 0.243705629816521, 0.243398357721994, 0.243029624980717, 0.242599429386111, 0.242107769165324, 0.24155464317143], [0.241407923874663, 0.241961050722937, 0.242452713414238, 0.242882912827661, 0.243251650322097, 0.243558927588949, 0.243804746468827, 0.243989108752087, 0.24411201598483, 0.244173469301379, 0.244173469301379, 0.24411201598483, 0.243989108752087, 0.243804746468827, 0.243558927588949, 0.243251650322097, 0.242882912827661, 0.242452713414238, 0.241961050722937, 0.241407923874663], [0.241261205902991, 0.241814333533832, 0.242305998487993, 0.24273620139926, 0.243104943247368, 0.243412225251851, 0.243658048739906, 0.243842415002574, 0.243965325154826, 0.244026780014718, 0.244026780014718, 0.243965325154826, 0.243842415002574, 0.243658048739906, 0.243412225251851, 0.243104943247368, 0.24273620139926, 0.242305998487993, 0.241814333533832, 0.241261205902991], [0.241114494986836, 0.241667623017863, 0.24215928912922, 0.2425894938293, 0.242958237904002, 0.24326552233158, 0.243511348176654, 0.24369571647484, 0.243818628120522, 0.24388008376991, 0.24388008376991, 0.243818628120522, 0.24369571647484, 0.243511348176654, 0.24326552233158, 0.242958237904002, 0.2425894938293, 0.24215928912922, 0.241667623017863, 0.241114494986836], [0.24096778407068, 0.241520912501894, 0.242012579770447, 0.24244278625934, 0.242811532560635, 0.243118819411309, 0.243364647613401, 0.243549017947107, 0.243671931086219, 0.243733387525101, 0.243733387525101, 0.243671931086219, 0.243549017947107, 0.243364647613401, 0.243118819411309, 0.242811532560635, 0.242442786259341, 0.242012579770447, 0.241520912501894, 0.24096778407068]]
temp_sol = [[2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98], [2.849766036993418, 2.849760178886377, 2.849754933371469, 2.849750313847679, 2.849746332114973, 2.84974299834416, 2.849740321050914, 2.849738307074019, 2.84973696155791, 2.849736287939525, 2.849736287939528, 2.849736961557921, 2.849738307074038, 2.849740321050939, 2.849742998344193, 2.849746332115012, 2.849750313847726, 2.849754933371523, 2.849760178886438, 2.849766036993486], [2.85002160238875, 2.850022143282418, 2.850022627613948, 2.850023054146176, 2.850023421789575, 2.85002372960504, 2.85002397680629, 2.850024162761873, 2.850024286996787, 2.850024349193681, 2.850024349193681, 2.850024286996787, 2.850024162761875, 2.850023976806292, 2.850023729605043, 2.850023421789578, 2.850023054146181, 2.850022627613954, 2.850022143282424, 2.850021602388757], [2.850007338781646, 2.850007522534495, 2.850007687072007, 2.85000783197389, 2.85000795687001, 2.850008061441333, 2.850008145420745, 2.850008208593727, 2.850008250798914, 2.850008271928492, 2.850008271928493, 2.850008250798913, 2.850008208593727, 2.850008145420745, 2.850008061441333, 2.85000795687001, 2.850007831973891, 2.850007687072008, 2.850007522534496, 2.850007338781647], [2.85000155642437, 2.850001595395063, 2.850001630290527, 2.850001661021623, 2.850001687509854, 2.850001709687559, 2.850001727498088, 2.850001740895946, 2.850001749846909, 2.850001754328114, 2.850001754328114, 2.850001749846909, 2.850001740895946, 2.850001727498087, 2.850001709687559, 2.850001687509854, 2.850001661021623, 2.850001630290526, 2.850001595395063, 2.850001556424369], [2.849999969542895, 2.849999968780291, 2.849999968097435, 2.849999967496069, 2.84999996697773, 2.849999966543743, 2.849999966195215, 2.849999965933037, 2.849999965757879, 2.849999965670187, 2.849999965670188, 2.84999996575788, 2.849999965933037, 2.849999966195215, 2.849999966543743, 2.84999996697773, 2.849999967496069, 2.849999968097435, 2.849999968780291, 2.849999969542895], [2.849999936330692, 2.849999934736502, 2.849999933309019, 2.849999932051889, 2.849999930968325, 2.849999930061092, 2.849999929332508, 2.849999928784437, 2.849999928418276, 2.849999928234961, 2.849999928234961, 2.849999928418276, 2.849999928784437, 2.849999929332508, 2.849999930061092, 2.849999930968324, 2.849999932051889, 2.849999933309019, 2.849999934736502, 2.849999936330693], [2.849999961457007, 2.849999960491944, 2.8499999596278, 2.849999958866781, 2.849999958210832, 2.849999957661627, 2.849999957220571, 2.849999956888789, 2.849999956667128, 2.849999956556156, 2.849999956556156, 2.849999956667128, 2.849999956888789, 2.84999995722057, 2.849999957661627, 2.849999958210831, 2.849999958866781, 2.849999959627799, 2.849999960491944, 2.849999961457007], [2.849999969643218, 2.849999968883127, 2.849999968202519, 2.849999967603135, 2.849999967086504, 2.849999966653946, 2.849999966306566, 2.849999966045252, 2.84999996587067, 2.849999965783268, 2.849999965783268, 2.84999996587067, 2.849999966045252, 2.849999966306566, 2.849999966653946, 2.849999967086503, 2.849999967603134, 2.849999968202519, 2.849999968883127, 2.849999969643219], [2.849999977829429, 2.849999977274309, 2.849999976777239, 2.849999976339488, 2.849999975962176, 2.849999975646264, 2.849999975392561, 2.849999975201715, 2.849999975074212, 2.849999975010379, 2.849999975010379, 2.849999975074212, 2.849999975201715, 2.849999975392561, 2.849999975646264, 2.849999975962175, 2.849999976339488, 2.849999976777239, 2.84999997727431, 2.84999997782943], [2.849999986015641, 2.849999985665493, 2.84999998535196, 2.849999985075843, 2.849999984837849, 2.849999984638584, 2.849999984478556, 2.849999984358178, 2.849999984277754, 2.849999984237491, 2.849999984237491, 2.849999984277754, 2.849999984358178, 2.849999984478556, 2.849999984638584, 2.849999984837848, 2.849999985075842, 2.849999985351959, 2.849999985665493, 2.849999986015641], [2.849999994201853, 2.849999994056676, 2.849999993926679, 2.849999993812196, 2.849999993713521, 2.849999993630902, 2.849999993564552, 2.849999993514642, 2.849999993481296, 2.849999993464602, 2.849999993464602, 2.849999993481296, 2.849999993514642, 2.849999993564552, 2.849999993630902, 2.84999999371352, 2.849999993812196, 2.849999993926679, 2.849999994056676, 2.849999994201853], [2.849999999555538, 2.849999999544409, 2.849999999534445, 2.849999999525669, 2.849999999518105, 2.849999999511772, 2.849999999506685, 2.84999999950286, 2.849999999500304, 2.849999999499023, 2.849999999499023, 2.849999999500304, 2.84999999950286, 2.849999999506685, 2.849999999511772, 2.849999999518105, 2.849999999525669, 2.849999999534444, 2.84999999954441, 2.849999999555538], [2.849999999689912, 2.849999999682148, 2.849999999675195, 2.849999999669073, 2.849999999663796, 2.849999999659377, 2.849999999655828, 2.84999999965316, 2.849999999651376, 2.849999999650484, 2.849999999650483, 2.849999999651376, 2.84999999965316, 2.849999999655828, 2.849999999659377, 2.849999999663796, 2.849999999669072, 2.849999999675195, 2.849999999682148, 2.849999999689912], [2.849999999824286, 2.849999999819885, 2.849999999815946, 2.849999999812477, 2.849999999809486, 2.849999999806983, 2.849999999804971, 2.84999999980346, 2.849999999802449, 2.849999999801943, 2.849999999801943, 2.849999999802449, 2.84999999980346, 2.849999999804971, 2.849999999806983, 2.849999999809486, 2.849999999812476, 2.849999999815946, 2.849999999819885, 2.849999999824286], [2.849999999958659, 2.849999999957623, 2.849999999956697, 2.84999999995588, 2.849999999955177, 2.849999999954588, 2.849999999954115, 2.849999999953759, 2.849999999953521, 2.849999999953403, 2.849999999953403, 2.849999999953521, 2.849999999953759, 2.849999999954115, 2.849999999954588, 2.849999999955177, 2.84999999995588, 2.849999999956697, 2.849999999957624, 2.849999999958659], [2.850000000093032, 2.850000000095362, 2.850000000097447, 2.850000000099284, 2.850000000100867, 2.850000000102193, 2.850000000103258, 2.850000000104059, 2.850000000104594, 2.850000000104862, 2.850000000104862, 2.850000000104594, 2.850000000104059, 2.850000000103258, 2.850000000102193, 2.850000000100867, 2.850000000099284, 2.850000000097447, 2.850000000095362, 2.850000000093032], [2.850000000152046, 2.850000000155853, 2.850000000159262, 2.850000000162264, 2.850000000164851, 2.850000000167018, 2.850000000168758, 2.850000000170066, 2.850000000170941, 2.850000000171379, 2.850000000171379, 2.850000000170941, 2.850000000170066, 2.850000000168758, 2.850000000167018, 2.850000000164852, 2.850000000162264, 2.850000000159262, 2.850000000155853, 2.850000000152046], [2.850000000119615, 2.85000000012261, 2.850000000125291, 2.850000000127653, 2.850000000129689, 2.850000000131394, 2.850000000132762, 2.850000000133791, 2.85000000013448, 2.850000000134824, 2.850000000134824, 2.85000000013448, 2.850000000133791, 2.850000000132762, 2.850000000131394, 2.850000000129689, 2.850000000127653, 2.850000000125291, 2.85000000012261, 2.850000000119615], [2.850000000087184, 2.850000000089366, 2.850000000091321, 2.850000000093043, 2.850000000094526, 2.850000000095769, 2.850000000096766, 2.850000000097516, 2.850000000098018, 2.850000000098269, 2.850000000098269, 2.850000000098019, 2.850000000097516, 2.850000000096766, 2.850000000095769, 2.850000000094526, 2.850000000093043, 2.850000000091321, 2.850000000089366, 2.850000000087184], [2.850000000054753, 2.850000000056123, 2.850000000057351, 2.850000000058432, 2.850000000059363, 2.850000000060144, 2.850000000060771, 2.850000000061241, 2.850000000061557, 2.850000000061714, 2.850000000061714, 2.850000000061557, 2.850000000061241, 2.850000000060771, 2.850000000060144, 2.850000000059364, 2.850000000058432, 2.850000000057351, 2.850000000056123, 2.850000000054753], [2.850000000022321, 2.85000000002288, 2.850000000023381, 2.850000000023821, 2.850000000024201, 2.850000000024519, 2.850000000024774, 2.850000000024966, 2.850000000025095, 2.85000000002516, 2.85000000002516, 2.850000000025095, 2.850000000024966, 2.850000000024774, 2.850000000024519, 2.850000000024201, 2.850000000023821, 2.850000000023381, 2.85000000002288, 2.850000000022321], [2.850000000006854, 2.850000000007025, 2.850000000007178, 2.850000000007313, 2.85000000000743, 2.850000000007527, 2.850000000007606, 2.850000000007666, 2.850000000007705, 2.850000000007725, 2.850000000007725, 2.850000000007705, 2.850000000007666, 2.850000000007606, 2.850000000007527, 2.85000000000743, 2.850000000007313, 2.850000000007179, 2.850000000007024, 2.850000000006854], [2.850000000005164, 2.850000000005293, 2.850000000005409, 2.85000000000551, 2.850000000005599, 2.850000000005672, 2.850000000005731, 2.850000000005776, 2.850000000005806, 2.850000000005821, 2.850000000005821, 2.850000000005806, 2.850000000005776, 2.850000000005731, 2.850000000005672, 2.850000000005599, 2.85000000000551, 2.850000000005409, 2.850000000005293, 2.850000000005164], [2.850000000003475, 2.850000000003562, 2.85000000000364, 2.850000000003708, 2.850000000003768, 2.850000000003817, 2.850000000003856, 2.850000000003887, 2.850000000003907, 2.850000000003917, 2.850000000003917, 2.850000000003907, 2.850000000003887, 2.850000000003856, 2.850000000003817, 2.850000000003767, 2.850000000003708, 2.85000000000364, 2.850000000003562, 2.850000000003475], [2.850000000001786, 2.85000000000183, 2.850000000001871, 2.850000000001905, 2.850000000001936, 2.850000000001961, 2.850000000001982, 2.850000000001998, 2.850000000002008, 2.850000000002013, 2.850000000002013, 2.850000000002008, 2.850000000001998, 2.850000000001982, 2.850000000001961, 2.850000000001936, 2.850000000001905, 2.850000000001871, 2.85000000000183, 2.850000000001786], [2.850000000000096, 2.850000000000099, 2.850000000000101, 2.850000000000102, 2.850000000000104, 2.850000000000106, 2.850000000000106, 2.850000000000108, 2.850000000000108, 2.850000000000108, 2.850000000000108, 2.850000000000108, 2.850000000000108, 2.850000000000106, 2.850000000000106, 2.850000000000104, 2.850000000000102, 2.850000000000101, 2.850000000000099, 2.850000000000096], [2.849999999999567, 2.849999999999556, 2.849999999999546, 2.849999999999537, 2.84999999999953, 2.849999999999524, 2.849999999999519, 2.849999999999516, 2.849999999999513, 2.849999999999512, 2.849999999999512, 2.849999999999513, 2.849999999999516, 2.849999999999519, 2.849999999999524, 2.84999999999953, 2.849999999999537, 2.849999999999546, 2.849999999999556, 2.849999999999568], [2.849999999999657, 2.849999999999648, 2.84999999999964, 2.849999999999633, 2.849999999999627, 2.849999999999623, 2.849999999999619, 2.849999999999616, 2.849999999999614, 2.849999999999613, 2.849999999999613, 2.849999999999614, 2.849999999999616, 2.849999999999619, 2.849999999999623, 2.849999999999627, 2.849999999999633, 2.849999999999641, 2.849999999999648, 2.849999999999657], [2.849999999999747, 2.84999999999974, 2.849999999999735, 2.849999999999729, 2.849999999999725, 2.849999999999722, 2.849999999999719, 2.849999999999716, 2.849999999999715, 2.849999999999715, 2.849999999999714, 2.849999999999715, 2.849999999999716, 2.849999999999719, 2.849999999999722, 2.849999999999725, 2.849999999999729, 2.849999999999735, 2.84999999999974, 2.849999999999747], [2.849999999999837, 2.849999999999832, 2.849999999999828, 2.849999999999825, 2.849999999999822, 2.84999999999982, 2.849999999999818, 2.849999999999817, 2.849999999999816, 2.849999999999816, 2.849999999999816, 2.849999999999816, 2.849999999999817, 2.849999999999818, 2.84999999999982, 2.849999999999822, 2.849999999999825, 2.849999999999828, 2.849999999999832, 2.849999999999837], [2.849999999999927, 2.849999999999925, 2.849999999999923, 2.849999999999921, 2.849999999999919, 2.849999999999919, 2.849999999999918, 2.849999999999917, 2.849999999999917, 2.849999999999917, 2.849999999999917, 2.849999999999917, 2.849999999999917, 2.849999999999918, 2.849999999999919, 2.849999999999919, 2.849999999999921, 2.849999999999923, 2.849999999999925, 2.849999999999927], [2.849999999999957, 2.849999999999956, 2.849999999999954, 2.849999999999953, 2.849999999999952, 2.849999999999952, 2.849999999999952, 2.849999999999951, 2.849999999999951, 2.849999999999951, 2.849999999999951, 2.849999999999951, 2.849999999999951, 2.849999999999952, 2.849999999999953, 2.849999999999953, 2.849999999999953, 2.849999999999954, 2.849999999999956, 2.849999999999957], [2.849999999999968, 2.849999999999966, 2.849999999999965, 2.849999999999965, 2.849999999999964, 2.849999999999964, 2.849999999999964, 2.849999999999963, 2.849999999999963, 2.849999999999963, 2.849999999999963, 2.849999999999963, 2.849999999999963, 2.849999999999964, 2.849999999999964, 2.849999999999964, 2.849999999999965, 2.849999999999966, 2.849999999999966, 2.849999999999968], [2.849999999999978, 2.849999999999977, 2.849999999999977, 2.849999999999976, 2.849999999999976, 2.849999999999976, 2.849999999999976, 2.849999999999975, 2.849999999999975, 2.849999999999975, 2.849999999999975, 2.849999999999975, 2.849999999999975, 2.849999999999976, 2.849999999999976, 2.849999999999976, 2.849999999999976, 2.849999999999977, 2.849999999999977, 2.849999999999978], [2.849999999999989, 2.849999999999988, 2.849999999999988, 2.849999999999988, 2.849999999999988, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999987, 2.849999999999988, 2.849999999999988, 2.849999999999988, 2.849999999999988, 2.849999999999988, 2.849999999999989], [2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999, 2.849999999999999], [2.850000000000001, 2.85, 2.850000000000001, 2.850000000000001, 2.85, 2.85, 2.85, 2.85, 2.850000000000001, 2.85, 2.85, 2.850000000000001, 2.85, 2.85, 2.850000000000001, 2.850000000000001, 2.85, 2.850000000000001, 2.850000000000001, 2.85], [2.85, 2.85, 2.85, 2.850000000000001, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.850000000000001, 2.85, 2.85, 2.850000000000001, 2.850000000000001, 2.85, 2.850000000000001, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85], [2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 2.85]]

####################################################

L = 2e-2 ##m thickness of wood
tend = 3600*360; ##s length of exposure time
node_x = 20;
node_t = 50;

m=0; ##0 is cartesian, 1 is cylindrical, 2 is spherical
##Heat Transfer - Inputs 
K = 0.2 ##W/m.K
rho = 600 ##kg/m^3
cp = 2500 ##J/kg-K
h = 5 ##convection heat transfer coefficient of air W/m^2.K
T_inf = 12+273##K

##Mass transfer- Inputs
Pa = 101325; ##air Pressure Pa
Cpa = 1010; ##Heat capacity of air J/kg.K
Ka = 0.0285; ##Thermal conductivity of air W/m.K
Ra = 287.058; ##Gas Constant J/kg.K
T = T_inf; ##K
Deff = 1e-10; ## moisture diffusivity, [m^2/s]
hm = 1E-10;

print(0.622*h*Deff**(2/3)/(Pa*Cpa**(1/3)*(Ka*Ra*T)**(2/3)))

Pwv_inf = 1170; ##pa
Pwv_sat = antoine_eqn(T-273);
RH = 0.6;

