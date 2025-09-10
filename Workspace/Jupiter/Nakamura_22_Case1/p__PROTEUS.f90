module p__PROTEUS
  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: get_number_of_reaction, get_number_of_all_species, get_number_of_var_species, get_number_of_fix_species, &
            get_all_species_name, get_var_species_name, get_fix_species_name, get_mass_of_species, &
            get_reaction_type_list, get_reactant_product_list, &
            get_number_of_wavelength_bins, get_wavelength_bin, get_solar_flux, &
            p__PROTEUS_source, p__PROTEUS_Jacobian 

contains


  integer function get_number_of_reaction()
    get_number_of_reaction = 218
    ! R1: H + hv -> H+ + e- 
    ! R2: H2 + hv -> H2+ + e- 
    ! R3: H2 + hv -> H+ + e- + H 
    ! R4: He + hv -> He+ + e- 
    ! R5: CH4 + hv -> CH4+ + e- 
    ! R6: CH4 + hv -> CH3+ + e- + H 
    ! R7: CH4 + hv -> H2+ + e- + products 
    ! R8: C2H2 + hv -> C2H2+ + e- 
    ! R9: C2H4 + hv -> C2H4+ + e- 
    ! R10: C2H4 + hv -> C2H3+ + e- + H 
    ! R11: C2H4 + hv -> C2H2+ + e- + products 
    ! R12: C2H4 + hv -> C2H+ + e- + products 
    ! R13: C2H6 + hv -> C2H6+ + e- 
    ! R14: C2H6 + hv -> C2H5+ + e- + H 
    ! R15: C2H6 + hv -> C2H4+ + e- + products 
    ! R16: C2H6 + hv -> C2H3+ + e- + products 
    ! R17: C2H6 + hv -> C2H2+ + e- + products 
    ! R18: H+ + e- -> H 
    ! R19: He+ + e- -> He 
    ! R20: HeH+ + e- -> H + He 
    ! R21: H2+ + e- -> H + H 
    ! R22: H3+ + e- -> H2 + H 
    ! R23: H3+ + e- -> H + H + H 
    ! R24: C+ + e- -> C 
    ! R25: CH+ + e- -> C + H 
    ! R26: CH2+ + e- -> CH + H 
    ! R27: CH3+ + e- -> CH2 + H 
    ! R28: CH4+ + e- -> CH3 + H 
    ! R29: CH4+ + e- -> CH2 + H + H 
    ! R30: CH5+ + e- -> CH2 + H + H2 
    ! R31: CH5+ + e- -> CH3 + H + H 
    ! R32: C2+ + e- -> C + C 
    ! R33: C2H+ + e- -> C2 + H 
    ! R34: C2H+ + e- -> CH + C 
    ! R35: C2H2+ + e- -> C2H + H 
    ! R36: C2H2+ + e- -> CH + CH 
    ! R37: C2H3+ + e- -> C2H2 + H 
    ! R38: C2H3+ + e- -> CH2 + CH 
    ! R39: C2H4+ + e- -> C2H3 + H 
    ! R40: C2H4+ + e- -> CH2 + CH2 
    ! R41: C2H5+ + e- -> C2H4 + H 
    ! R42: C2H5+ + e- -> CH3 + CH2 
    ! R43: C2H6+ + e- -> C2H5 + H 
    ! R44: C2H6+ + e- -> CH3 + CH3 
    ! R45: C2H7+ + e- -> C2H6 + H 
    ! R46: C3H+ + e- -> products 
    ! R47: C3H2+ + e- -> products 
    ! R48: C3H3+ + e- -> products 
    ! R49: C3H4+ + e- -> products 
    ! R50: C3H5+ + e- -> products 
    ! R51: C3H6+ + e- -> products 
    ! R52: C3H7+ + e- -> products 
    ! R53: C3H8+ + e- -> products 
    ! R54: C3H9+ + e- -> products 
    ! R55: C4H+ + e- -> products 
    ! R56: C4H2+ + e- -> products 
    ! R57: C4H3+ + e- -> products 
    ! R58: C4H5+ + e- -> products 
    ! R59: C4H7+ + e- -> products 
    ! R60: C4H9+ + e- -> products 
    ! R61: H2+ + H2 -> H3+ + H 
    ! R62: H2+ + H -> H+ + H2 
    ! R63: H2+ + He -> HeH+ + H 
    ! R64: H2+ + CH4 -> CH5+ + H 
    ! R65: H2+ + CH4 -> CH4+ + H2 
    ! R66: H2+ + CH4 -> CH3+ + H + H2 
    ! R67: H2+ + C2H2 -> C2H3+ + H 
    ! R68: H2+ + C2H2 -> C2H2+ + H2 
    ! R69: H2+ + C2H6 -> C2H6+ + H2 
    ! R70: H2+ + C2H6 -> C2H5+ + H + H2 
    ! R71: H2+ + C2H6 -> C2H4+ + H2 + H2 
    ! R72: H2+ + C2H6 -> C2H3+ + H + H2 + H2 
    ! R73: H2+ + C2H6 -> C2H2+ + H2 + H2 + H2 
    ! R74: He+ + H2 -> H2+ + He 
    ! R75: He+ + H2 -> H+ + H + He 
    ! R76: He+ + H2(v2) -> H+ + H + He 
    ! R77: He+ + CH4 -> H+ + CH3 + He 
    ! R78: He+ + CH4 -> CH+ + H2 + H + He 
    ! R79: He+ + CH4 -> CH2+ + H2 + He 
    ! R80: He+ + CH4 -> CH3+ + H + He 
    ! R81: He+ + CH4 -> CH4+ + He 
    ! R82: He+ + C2H2 -> C2H2+ + He 
    ! R83: He+ + C2H2 -> C2H+ + H + He 
    ! R84: He+ + C2H2 -> C2+ + H2 + He 
    ! R85: He+ + C2H2 -> CH+ + CH + He 
    ! R86: He+ + C2H4 -> C2H4+ + He 
    ! R87: He+ + C2H4 -> C2H3+ + H + He 
    ! R88: He+ + C2H4 -> C2H2+ + H2 + He 
    ! R89: He+ + C2H4 -> C2H+ + H + H2 + He 
    ! R90: He+ + C2H4 -> CH2+ + CH2 + He 
    ! R91: He+ + C2H6 -> C2H4+ + H2 + He 
    ! R92: He+ + C2H6 -> C2H3+ + H + H2 + He 
    ! R93: He+ + C2H6 -> C2H2+ + H2 + H2 + He 
    ! R94: H+ + CH4 -> CH3+ + H2 
    ! R95: H+ + CH4 -> CH4+ + H 
    ! R96: H+ + C2H6 -> C2H3+ + H2 + H2 
    ! R97: H+ + C2H6 -> C2H4+ + H + H2 
    ! R98: H+ + C2H6 -> C2H5+ + H2 
    ! R99: H+ + H2(v4) -> H2+ + H 
    ! R100: HeH+ + H2 -> H3+ + He 
    ! R101: HeH+ + H -> H2+ + He 
    ! R102: HeH+ + C2H4 -> C2H4+ + H + He 
    ! R103: HeH+ + C2H4 -> C2H3+ + H2 + He 
    ! R104: HeH+ + C2H6 -> C2H3+ + H2 + H2 + He 
    ! R105: HeH+ + C2H6 -> C2H5+ + H2 + He 
    ! R106: H3+ + CH4 -> CH5+ + H2 
    ! R107: H3+ + C2H2 -> C2H3+ + H2 
    ! R108: H3+ + C2H4 -> C2H5+ + H2 
    ! R109: H3+ + C2H4 -> C2H3+ + H2 + H2 
    ! R110: H3+ + C2H6 -> C2H5+ + H2 + H2 
    ! R111: C+ + CH4 -> C2H2+ + H2 
    ! R112: C+ + CH4 -> C2H3+ + H 
    ! R113: C+ + C2H2 -> C3H+ + H 
    ! R114: C+ + C2H4 -> C2H3+ + CH 
    ! R115: C+ + C2H4 -> C2H4+ + C 
    ! R116: C+ + C2H4 -> C3H+ + H + H2 
    ! R117: C+ + C2H4 -> C3H2+ + H2 
    ! R118: C+ + C2H4 -> C3H3+ + H 
    ! R119: C+ + C2H6 -> C2H2+ + CH4 
    ! R120: C+ + C2H6 -> C2H3+ + CH3 
    ! R121: C+ + C2H6 -> C2H5+ + CH 
    ! R122: C+ + C2H6 -> C3H3+ + H + H2 
    ! R123: CH+ + H -> C+ + H2 
    ! R124: CH+ + H2 -> CH2+ + H 
    ! R125: CH+ + CH4 -> C2H4+ + H 
    ! R126: CH+ + CH4 -> C2H3+ + H2 
    ! R127: CH+ + CH4 -> C2H2+ + H + H2 
    ! R128: CH+ + C2H2 -> C3H2+ + H 
    ! R129: CH2+ + H2 -> CH3+ + H 
    ! R130: CH2+ + CH4 -> C2H5+ + H 
    ! R131: CH2+ + CH4 -> C2H4+ + H2 
    ! R132: CH2+ + C2H2 -> C3H3+ + H 
    ! R133: CH3+ + CH4 -> C2H5+ + H2 
    ! R134: CH3+ + C2H2 -> C3H3+ + H2 
    ! R135: CH3+ + C2H4 -> C2H3+ + CH4 
    ! R136: CH3+ + C2H4 -> C3H3+ + H2 + H2 
    ! R137: CH3+ + C2H4 -> C3H5+ + H2 
    ! R138: CH3+ + C2H6 -> C2H5+ + CH4 
    ! R139: CH3+ + C2H6 -> C3H5+ + H2 + H2 
    ! R140: CH3+ + C2H6 -> C3H7+ + H2 
    ! R141: CH4+ + H2 -> CH5+ + H 
    ! R142: CH4+ + CH4 -> CH5+ + CH3 
    ! R143: CH4+ + C2H2 -> C2H2+ + CH4 
    ! R144: CH4+ + C2H2 -> C2H3+ + CH3 
    ! R145: CH4+ + C2H2 -> C3H3+ + H + H2 
    ! R146: CH4+ + C2H4 -> C2H4+ + CH4 
    ! R147: CH4+ + C2H4 -> C2H5+ + CH3 
    ! R148: CH4+ + C2H4 -> C3H5+ + H + H2 
    ! R149: CH4+ + C2H6 -> C2H4+ + CH4 + H2 
    ! R150: CH5+ + H -> CH4+ + H2 
    ! R151: CH5+ + C2H2 -> C2H3+ + CH4 
    ! R152: CH5+ + C2H4 -> C2H5+ + CH4 
    ! R153: CH5+ + C2H6 -> C2H5+ + CH4 + H2 
    ! R154: CH5+ + C2H6 -> C2H7+ + CH4 
    ! R155: C2+ + H2 -> C2H+ + H 
    ! R156: C2+ + CH4 -> C2H+ + CH3 
    ! R157: C2+ + CH4 -> C2H2+ + CH2 
    ! R158: C2+ + CH4 -> C3H+ + H + H2 
    ! R159: C2+ + CH4 -> C3H2+ + H2 
    ! R160: C2+ + CH4 -> C3H3+ + H 
    ! R161: C2+ + C2H2 -> C4H+ + H 
    ! R162: C2+ + C2H4 -> products 
    ! R163: C2H+ + H2 -> C2H2+ + H 
    ! R164: C2H+ + CH4 -> C2H2+ + CH3 
    ! R165: C2H+ + CH4 -> C3H3+ + H2 
    ! R166: C2H+ + CH4 -> C3H4+ + H 
    ! R167: C2H+ + CH4 -> C3H5+ 
    ! R168: C2H+ + C2H2 -> C4H2+ + H 
    ! R169: C2H+ + C2H4 -> products 
    ! R170: C2H2+ + H2 -> C2H3+ + H 
    ! R171: C2H2+ + CH4 -> C3H4+ + H2 
    ! R172: C2H2+ + CH4 -> C3H5+ + H 
    ! R173: C2H2+ + C2H2 -> C4H2+ + H2 
    ! R174: C2H2+ + C2H2 -> C4H3+ + H 
    ! R175: C2H2+ + C2H4 -> C2H4+ + C2H2 
    ! R176: C2H2+ + C2H4 -> C3H3+ + CH3 
    ! R177: C2H2+ + C2H4 -> C4H5+ + H 
    ! R178: C2H2+ + C2H6 -> C2H4+ + C2H4 
    ! R179: C2H2+ + C2H6 -> C2H5+ + C2H3 
    ! R180: C2H2+ + C2H6 -> C3H3+ + CH3 + H2 
    ! R181: C2H2+ + C2H6 -> C3H4+ + CH4 
    ! R182: C2H2+ + C2H6 -> C3H5+ + CH3 
    ! R183: C2H2+ + C2H6 -> C4H5+ + H + H2 
    ! R184: C2H2+ + C2H6 -> C4H7+ + H 
    ! R185: C2H3+ + H -> C2H2+ + H2 
    ! R186: C2H3+ + CH4 -> C3H5+ + H2 
    ! R187: C2H3+ + C2H2 -> C4H3+ + H2 
    ! R188: C2H3+ + C2H4 -> C2H5+ + C2H2 
    ! R189: C2H3+ + C2H6 -> C2H5+ + C2H4 
    ! R190: C2H3+ + C2H6 -> C3H5+ + CH4 
    ! R191: C2H3+ + C2H6 -> C4H7+ + H2 
    ! R192: C2H4+ + H -> C2H3+ + H2 
    ! R193: C2H4+ + C2H2 -> C3H3+ + CH3 
    ! R194: C2H4+ + C2H2 -> C4H5+ + H 
    ! R195: C2H4+ + C2H4 -> C3H5+ + CH3 
    ! R196: C2H4+ + C2H4 -> C4H7+ + H 
    ! R197: C2H4+ + C2H6 -> C3H6+ + CH4 
    ! R198: C2H4+ + C2H6 -> C3H7+ + CH3 
    ! R199: C2H5+ + H -> C2H4+ + H2 
    ! R200: C2H5+ + CH4 -> C3H7+ + H2 
    ! R201: C2H5+ + C2H2 -> C3H3+ + CH4 
    ! R202: C2H5+ + C2H2 -> C4H5+ + H2 
    ! R203: C2H5+ + C2H4 -> C3H5+ + CH4 
    ! R204: C2H5+ + C2H6 -> C4H9+ + H2 
    ! R205: C2H6+ + H -> C2H5+ + H2 
    ! R206: C2H6+ + C2H2 -> C2H5+ + C2H3 
    ! R207: C2H6+ + C2H2 -> C3H5+ + CH3 
    ! R208: C2H6+ + C2H2 -> C4H7+ + H 
    ! R209: C2H6+ + C2H4 -> C2H4+ + C2H6 
    ! R210: C2H6+ + C2H6 -> C3H8+ + CH4 
    ! R211: C2H6+ + C2H6 -> C3H9+ + CH3 
    ! R212: C2H7+ + C2H2 -> C2H3+ + C2H6 
    ! R213: C2H7+ + C2H4 -> C2H5+ + C2H6 
    ! R214: H+ + H2 + H2 -> H3+ + H2 
    ! R215: CH3+ + H2 + H2 -> CH5+ + H2 
    ! R216: CH3+ + H2 + He -> CH5+ + He 
    ! R217: C2H2+ + H2 + He -> C2H4+ + He 
    ! R218: H + H + H2 -> H2 + H2 
  end function


  integer function get_number_of_all_species()
    get_number_of_all_species = 52
  end function


  integer function get_number_of_var_species()
    get_number_of_var_species = 35
  end function


  integer function get_number_of_fix_species()
    get_number_of_fix_species = 17
  end function


  function get_all_species_name()
    implicit none
    character(len=256), dimension(52) :: get_all_species_name
    get_all_species_name(1) = 'H'
    get_all_species_name(2) = 'H+'
    get_all_species_name(3) = 'e-'
    get_all_species_name(4) = 'H2'
    get_all_species_name(5) = 'H2+'
    get_all_species_name(6) = 'He'
    get_all_species_name(7) = 'He+'
    get_all_species_name(8) = 'CH4'
    get_all_species_name(9) = 'CH4+'
    get_all_species_name(10) = 'CH3+'
    get_all_species_name(11) = 'C2H2'
    get_all_species_name(12) = 'C2H2+'
    get_all_species_name(13) = 'C2H4'
    get_all_species_name(14) = 'C2H4+'
    get_all_species_name(15) = 'C2H3+'
    get_all_species_name(16) = 'C2H+'
    get_all_species_name(17) = 'C2H6'
    get_all_species_name(18) = 'C2H6+'
    get_all_species_name(19) = 'C2H5+'
    get_all_species_name(20) = 'HeH+'
    get_all_species_name(21) = 'H3+'
    get_all_species_name(22) = 'C+'
    get_all_species_name(23) = 'C'
    get_all_species_name(24) = 'CH+'
    get_all_species_name(25) = 'CH2+'
    get_all_species_name(26) = 'CH'
    get_all_species_name(27) = 'CH2'
    get_all_species_name(28) = 'CH3'
    get_all_species_name(29) = 'CH5+'
    get_all_species_name(30) = 'C2+'
    get_all_species_name(31) = 'C2'
    get_all_species_name(32) = 'C2H'
    get_all_species_name(33) = 'C2H3'
    get_all_species_name(34) = 'C2H5'
    get_all_species_name(35) = 'C2H7+'
    get_all_species_name(36) = 'C3H+'
    get_all_species_name(37) = 'C3H2+'
    get_all_species_name(38) = 'C3H3+'
    get_all_species_name(39) = 'C3H4+'
    get_all_species_name(40) = 'C3H5+'
    get_all_species_name(41) = 'C3H6+'
    get_all_species_name(42) = 'C3H7+'
    get_all_species_name(43) = 'C3H8+'
    get_all_species_name(44) = 'C3H9+'
    get_all_species_name(45) = 'C4H+'
    get_all_species_name(46) = 'C4H2+'
    get_all_species_name(47) = 'C4H3+'
    get_all_species_name(48) = 'C4H5+'
    get_all_species_name(49) = 'C4H7+'
    get_all_species_name(50) = 'C4H9+'
    get_all_species_name(51) = 'H2(v2)'
    get_all_species_name(52) = 'H2(v4)'
  end function


  function get_var_species_name()
    implicit none
    character(len=256), dimension(35) :: get_var_species_name
    get_var_species_name(1) = 'H'
    get_var_species_name(2) = 'H+'
    get_var_species_name(3) = 'H2+'
    get_var_species_name(4) = 'He+'
    get_var_species_name(5) = 'CH4+'
    get_var_species_name(6) = 'CH3+'
    get_var_species_name(7) = 'C2H2+'
    get_var_species_name(8) = 'C2H4+'
    get_var_species_name(9) = 'C2H3+'
    get_var_species_name(10) = 'C2H+'
    get_var_species_name(11) = 'C2H6+'
    get_var_species_name(12) = 'C2H5+'
    get_var_species_name(13) = 'HeH+'
    get_var_species_name(14) = 'H3+'
    get_var_species_name(15) = 'C+'
    get_var_species_name(16) = 'CH+'
    get_var_species_name(17) = 'CH2+'
    get_var_species_name(18) = 'CH5+'
    get_var_species_name(19) = 'C2+'
    get_var_species_name(20) = 'C2H7+'
    get_var_species_name(21) = 'C3H+'
    get_var_species_name(22) = 'C3H2+'
    get_var_species_name(23) = 'C3H3+'
    get_var_species_name(24) = 'C3H4+'
    get_var_species_name(25) = 'C3H5+'
    get_var_species_name(26) = 'C3H6+'
    get_var_species_name(27) = 'C3H7+'
    get_var_species_name(28) = 'C3H8+'
    get_var_species_name(29) = 'C3H9+'
    get_var_species_name(30) = 'C4H+'
    get_var_species_name(31) = 'C4H2+'
    get_var_species_name(32) = 'C4H3+'
    get_var_species_name(33) = 'C4H5+'
    get_var_species_name(34) = 'C4H7+'
    get_var_species_name(35) = 'C4H9+'
  end function


  function get_fix_species_name()
    implicit none
    character(len=256), dimension(17) :: get_fix_species_name
    get_fix_species_name(1) = 'e-'
    get_fix_species_name(2) = 'H2'
    get_fix_species_name(3) = 'He'
    get_fix_species_name(4) = 'CH4'
    get_fix_species_name(5) = 'C2H2'
    get_fix_species_name(6) = 'C2H4'
    get_fix_species_name(7) = 'C2H6'
    get_fix_species_name(8) = 'C'
    get_fix_species_name(9) = 'CH'
    get_fix_species_name(10) = 'CH2'
    get_fix_species_name(11) = 'CH3'
    get_fix_species_name(12) = 'C2'
    get_fix_species_name(13) = 'C2H'
    get_fix_species_name(14) = 'C2H3'
    get_fix_species_name(15) = 'C2H5'
    get_fix_species_name(16) = 'H2(v2)'
    get_fix_species_name(17) = 'H2(v4)'
  end function


  function get_mass_of_species()
    implicit none
    integer, parameter :: nsp = 52
    real(dp), dimension(1:nsp) :: get_mass_of_species
    real(dp), parameter :: m_u = 1.660538921e-27_dp ! kg
    get_mass_of_species = 0.0_dp
    get_mass_of_species(1) = 1.0_dp * m_u ! H
    get_mass_of_species(2) = 1.0_dp * m_u ! H+
    get_mass_of_species(3) = 0.00054858_dp * m_u ! e-
    get_mass_of_species(4) = 2.0_dp * m_u ! H2
    get_mass_of_species(5) = 2.0_dp * m_u ! H2+
    get_mass_of_species(6) = 4.0_dp * m_u ! He
    get_mass_of_species(7) = 4.0_dp * m_u ! He+
    get_mass_of_species(8) = 16.0_dp * m_u ! CH4
    get_mass_of_species(9) = 16.0_dp * m_u ! CH4+
    get_mass_of_species(10) = 15.0_dp * m_u ! CH3+
    get_mass_of_species(11) = 26.0_dp * m_u ! C2H2
    get_mass_of_species(12) = 26.0_dp * m_u ! C2H2+
    get_mass_of_species(13) = 28.0_dp * m_u ! C2H4
    get_mass_of_species(14) = 28.0_dp * m_u ! C2H4+
    get_mass_of_species(15) = 27.0_dp * m_u ! C2H3+
    get_mass_of_species(16) = 25.0_dp * m_u ! C2H+
    get_mass_of_species(17) = 30.0_dp * m_u ! C2H6
    get_mass_of_species(18) = 30.0_dp * m_u ! C2H6+
    get_mass_of_species(19) = 29.0_dp * m_u ! C2H5+
    get_mass_of_species(20) = 5.0_dp * m_u ! HeH+
    get_mass_of_species(21) = 3.0_dp * m_u ! H3+
    get_mass_of_species(22) = 12.0_dp * m_u ! C+
    get_mass_of_species(23) = 12.0_dp * m_u ! C
    get_mass_of_species(24) = 13.0_dp * m_u ! CH+
    get_mass_of_species(25) = 14.0_dp * m_u ! CH2+
    get_mass_of_species(26) = 13.0_dp * m_u ! CH
    get_mass_of_species(27) = 14.0_dp * m_u ! CH2
    get_mass_of_species(28) = 15.0_dp * m_u ! CH3
    get_mass_of_species(29) = 17.0_dp * m_u ! CH5+
    get_mass_of_species(30) = 24.0_dp * m_u ! C2+
    get_mass_of_species(31) = 24.0_dp * m_u ! C2
    get_mass_of_species(32) = 25.0_dp * m_u ! C2H
    get_mass_of_species(33) = 27.0_dp * m_u ! C2H3
    get_mass_of_species(34) = 29.0_dp * m_u ! C2H5
    get_mass_of_species(35) = 31.0_dp * m_u ! C2H7+
    get_mass_of_species(36) = 37.0_dp * m_u ! C3H+
    get_mass_of_species(37) = 38.0_dp * m_u ! C3H2+
    get_mass_of_species(38) = 39.0_dp * m_u ! C3H3+
    get_mass_of_species(39) = 40.0_dp * m_u ! C3H4+
    get_mass_of_species(40) = 41.0_dp * m_u ! C3H5+
    get_mass_of_species(41) = 42.0_dp * m_u ! C3H6+
    get_mass_of_species(42) = 43.0_dp * m_u ! C3H7+
    get_mass_of_species(43) = 44.0_dp * m_u ! C3H8+
    get_mass_of_species(44) = 45.0_dp * m_u ! C3H9+
    get_mass_of_species(45) = 49.0_dp * m_u ! C4H+
    get_mass_of_species(46) = 50.0_dp * m_u ! C4H2+
    get_mass_of_species(47) = 51.0_dp * m_u ! C4H3+
    get_mass_of_species(48) = 53.0_dp * m_u ! C4H5+
    get_mass_of_species(49) = 55.0_dp * m_u ! C4H7+
    get_mass_of_species(50) = 57.0_dp * m_u ! C4H9+
    get_mass_of_species(51) = 2.0_dp * m_u ! H2(v2)
    get_mass_of_species(52) = 2.0_dp * m_u ! H2(v4)
  end function get_mass_of_species


  subroutine get_reactant_product_list(reactant_list, product_list)
    implicit none
    integer, parameter :: nch = 218
    integer, intent(out) :: reactant_list(1:nch,0:20)
    integer, intent(out) :: product_list(1:nch,0:20)

    ! reactant list ----------------------------
    reactant_list(1,0:1) = (/1,1/)
    reactant_list(2,0:1) = (/1,4/)
    reactant_list(3,0:1) = (/1,4/)
    reactant_list(4,0:1) = (/1,6/)
    reactant_list(5,0:1) = (/1,8/)
    reactant_list(6,0:1) = (/1,8/)
    reactant_list(7,0:1) = (/1,8/)
    reactant_list(8,0:1) = (/1,11/)
    reactant_list(9,0:1) = (/1,13/)
    reactant_list(10,0:1) = (/1,13/)
    reactant_list(11,0:1) = (/1,13/)
    reactant_list(12,0:1) = (/1,13/)
    reactant_list(13,0:1) = (/1,17/)
    reactant_list(14,0:1) = (/1,17/)
    reactant_list(15,0:1) = (/1,17/)
    reactant_list(16,0:1) = (/1,17/)
    reactant_list(17,0:1) = (/1,17/)
    reactant_list(18,0:2) = (/2,2,3/)
    reactant_list(19,0:2) = (/2,7,3/)
    reactant_list(20,0:2) = (/2,20,3/)
    reactant_list(21,0:2) = (/2,5,3/)
    reactant_list(22,0:2) = (/2,21,3/)
    reactant_list(23,0:2) = (/2,21,3/)
    reactant_list(24,0:2) = (/2,22,3/)
    reactant_list(25,0:2) = (/2,24,3/)
    reactant_list(26,0:2) = (/2,25,3/)
    reactant_list(27,0:2) = (/2,10,3/)
    reactant_list(28,0:2) = (/2,9,3/)
    reactant_list(29,0:2) = (/2,9,3/)
    reactant_list(30,0:2) = (/2,29,3/)
    reactant_list(31,0:2) = (/2,29,3/)
    reactant_list(32,0:2) = (/2,30,3/)
    reactant_list(33,0:2) = (/2,16,3/)
    reactant_list(34,0:2) = (/2,16,3/)
    reactant_list(35,0:2) = (/2,12,3/)
    reactant_list(36,0:2) = (/2,12,3/)
    reactant_list(37,0:2) = (/2,15,3/)
    reactant_list(38,0:2) = (/2,15,3/)
    reactant_list(39,0:2) = (/2,14,3/)
    reactant_list(40,0:2) = (/2,14,3/)
    reactant_list(41,0:2) = (/2,19,3/)
    reactant_list(42,0:2) = (/2,19,3/)
    reactant_list(43,0:2) = (/2,18,3/)
    reactant_list(44,0:2) = (/2,18,3/)
    reactant_list(45,0:2) = (/2,35,3/)
    reactant_list(46,0:2) = (/2,36,3/)
    reactant_list(47,0:2) = (/2,37,3/)
    reactant_list(48,0:2) = (/2,38,3/)
    reactant_list(49,0:2) = (/2,39,3/)
    reactant_list(50,0:2) = (/2,40,3/)
    reactant_list(51,0:2) = (/2,41,3/)
    reactant_list(52,0:2) = (/2,42,3/)
    reactant_list(53,0:2) = (/2,43,3/)
    reactant_list(54,0:2) = (/2,44,3/)
    reactant_list(55,0:2) = (/2,45,3/)
    reactant_list(56,0:2) = (/2,46,3/)
    reactant_list(57,0:2) = (/2,47,3/)
    reactant_list(58,0:2) = (/2,48,3/)
    reactant_list(59,0:2) = (/2,49,3/)
    reactant_list(60,0:2) = (/2,50,3/)
    reactant_list(61,0:2) = (/2,5,4/)
    reactant_list(62,0:2) = (/2,5,1/)
    reactant_list(63,0:2) = (/2,5,6/)
    reactant_list(64,0:2) = (/2,5,8/)
    reactant_list(65,0:2) = (/2,5,8/)
    reactant_list(66,0:2) = (/2,5,8/)
    reactant_list(67,0:2) = (/2,5,11/)
    reactant_list(68,0:2) = (/2,5,11/)
    reactant_list(69,0:2) = (/2,5,17/)
    reactant_list(70,0:2) = (/2,5,17/)
    reactant_list(71,0:2) = (/2,5,17/)
    reactant_list(72,0:2) = (/2,5,17/)
    reactant_list(73,0:2) = (/2,5,17/)
    reactant_list(74,0:2) = (/2,7,4/)
    reactant_list(75,0:2) = (/2,7,4/)
    reactant_list(76,0:2) = (/2,7,51/)
    reactant_list(77,0:2) = (/2,7,8/)
    reactant_list(78,0:2) = (/2,7,8/)
    reactant_list(79,0:2) = (/2,7,8/)
    reactant_list(80,0:2) = (/2,7,8/)
    reactant_list(81,0:2) = (/2,7,8/)
    reactant_list(82,0:2) = (/2,7,11/)
    reactant_list(83,0:2) = (/2,7,11/)
    reactant_list(84,0:2) = (/2,7,11/)
    reactant_list(85,0:2) = (/2,7,11/)
    reactant_list(86,0:2) = (/2,7,13/)
    reactant_list(87,0:2) = (/2,7,13/)
    reactant_list(88,0:2) = (/2,7,13/)
    reactant_list(89,0:2) = (/2,7,13/)
    reactant_list(90,0:2) = (/2,7,13/)
    reactant_list(91,0:2) = (/2,7,17/)
    reactant_list(92,0:2) = (/2,7,17/)
    reactant_list(93,0:2) = (/2,7,17/)
    reactant_list(94,0:2) = (/2,2,8/)
    reactant_list(95,0:2) = (/2,2,8/)
    reactant_list(96,0:2) = (/2,2,17/)
    reactant_list(97,0:2) = (/2,2,17/)
    reactant_list(98,0:2) = (/2,2,17/)
    reactant_list(99,0:2) = (/2,2,52/)
    reactant_list(100,0:2) = (/2,20,4/)
    reactant_list(101,0:2) = (/2,20,1/)
    reactant_list(102,0:2) = (/2,20,13/)
    reactant_list(103,0:2) = (/2,20,13/)
    reactant_list(104,0:2) = (/2,20,17/)
    reactant_list(105,0:2) = (/2,20,17/)
    reactant_list(106,0:2) = (/2,21,8/)
    reactant_list(107,0:2) = (/2,21,11/)
    reactant_list(108,0:2) = (/2,21,13/)
    reactant_list(109,0:2) = (/2,21,13/)
    reactant_list(110,0:2) = (/2,21,17/)
    reactant_list(111,0:2) = (/2,22,8/)
    reactant_list(112,0:2) = (/2,22,8/)
    reactant_list(113,0:2) = (/2,22,11/)
    reactant_list(114,0:2) = (/2,22,13/)
    reactant_list(115,0:2) = (/2,22,13/)
    reactant_list(116,0:2) = (/2,22,13/)
    reactant_list(117,0:2) = (/2,22,13/)
    reactant_list(118,0:2) = (/2,22,13/)
    reactant_list(119,0:2) = (/2,22,17/)
    reactant_list(120,0:2) = (/2,22,17/)
    reactant_list(121,0:2) = (/2,22,17/)
    reactant_list(122,0:2) = (/2,22,17/)
    reactant_list(123,0:2) = (/2,24,1/)
    reactant_list(124,0:2) = (/2,24,4/)
    reactant_list(125,0:2) = (/2,24,8/)
    reactant_list(126,0:2) = (/2,24,8/)
    reactant_list(127,0:2) = (/2,24,8/)
    reactant_list(128,0:2) = (/2,24,11/)
    reactant_list(129,0:2) = (/2,25,4/)
    reactant_list(130,0:2) = (/2,25,8/)
    reactant_list(131,0:2) = (/2,25,8/)
    reactant_list(132,0:2) = (/2,25,11/)
    reactant_list(133,0:2) = (/2,10,8/)
    reactant_list(134,0:2) = (/2,10,11/)
    reactant_list(135,0:2) = (/2,10,13/)
    reactant_list(136,0:2) = (/2,10,13/)
    reactant_list(137,0:2) = (/2,10,13/)
    reactant_list(138,0:2) = (/2,10,17/)
    reactant_list(139,0:2) = (/2,10,17/)
    reactant_list(140,0:2) = (/2,10,17/)
    reactant_list(141,0:2) = (/2,9,4/)
    reactant_list(142,0:2) = (/2,9,8/)
    reactant_list(143,0:2) = (/2,9,11/)
    reactant_list(144,0:2) = (/2,9,11/)
    reactant_list(145,0:2) = (/2,9,11/)
    reactant_list(146,0:2) = (/2,9,13/)
    reactant_list(147,0:2) = (/2,9,13/)
    reactant_list(148,0:2) = (/2,9,13/)
    reactant_list(149,0:2) = (/2,9,17/)
    reactant_list(150,0:2) = (/2,29,1/)
    reactant_list(151,0:2) = (/2,29,11/)
    reactant_list(152,0:2) = (/2,29,13/)
    reactant_list(153,0:2) = (/2,29,17/)
    reactant_list(154,0:2) = (/2,29,17/)
    reactant_list(155,0:2) = (/2,30,4/)
    reactant_list(156,0:2) = (/2,30,8/)
    reactant_list(157,0:2) = (/2,30,8/)
    reactant_list(158,0:2) = (/2,30,8/)
    reactant_list(159,0:2) = (/2,30,8/)
    reactant_list(160,0:2) = (/2,30,8/)
    reactant_list(161,0:2) = (/2,30,11/)
    reactant_list(162,0:2) = (/2,30,13/)
    reactant_list(163,0:2) = (/2,16,4/)
    reactant_list(164,0:2) = (/2,16,8/)
    reactant_list(165,0:2) = (/2,16,8/)
    reactant_list(166,0:2) = (/2,16,8/)
    reactant_list(167,0:2) = (/2,16,8/)
    reactant_list(168,0:2) = (/2,16,11/)
    reactant_list(169,0:2) = (/2,16,13/)
    reactant_list(170,0:2) = (/2,12,4/)
    reactant_list(171,0:2) = (/2,12,8/)
    reactant_list(172,0:2) = (/2,12,8/)
    reactant_list(173,0:2) = (/2,12,11/)
    reactant_list(174,0:2) = (/2,12,11/)
    reactant_list(175,0:2) = (/2,12,13/)
    reactant_list(176,0:2) = (/2,12,13/)
    reactant_list(177,0:2) = (/2,12,13/)
    reactant_list(178,0:2) = (/2,12,17/)
    reactant_list(179,0:2) = (/2,12,17/)
    reactant_list(180,0:2) = (/2,12,17/)
    reactant_list(181,0:2) = (/2,12,17/)
    reactant_list(182,0:2) = (/2,12,17/)
    reactant_list(183,0:2) = (/2,12,17/)
    reactant_list(184,0:2) = (/2,12,17/)
    reactant_list(185,0:2) = (/2,15,1/)
    reactant_list(186,0:2) = (/2,15,8/)
    reactant_list(187,0:2) = (/2,15,11/)
    reactant_list(188,0:2) = (/2,15,13/)
    reactant_list(189,0:2) = (/2,15,17/)
    reactant_list(190,0:2) = (/2,15,17/)
    reactant_list(191,0:2) = (/2,15,17/)
    reactant_list(192,0:2) = (/2,14,1/)
    reactant_list(193,0:2) = (/2,14,11/)
    reactant_list(194,0:2) = (/2,14,11/)
    reactant_list(195,0:2) = (/2,14,13/)
    reactant_list(196,0:2) = (/2,14,13/)
    reactant_list(197,0:2) = (/2,14,17/)
    reactant_list(198,0:2) = (/2,14,17/)
    reactant_list(199,0:2) = (/2,19,1/)
    reactant_list(200,0:2) = (/2,19,8/)
    reactant_list(201,0:2) = (/2,19,11/)
    reactant_list(202,0:2) = (/2,19,11/)
    reactant_list(203,0:2) = (/2,19,13/)
    reactant_list(204,0:2) = (/2,19,17/)
    reactant_list(205,0:2) = (/2,18,1/)
    reactant_list(206,0:2) = (/2,18,11/)
    reactant_list(207,0:2) = (/2,18,11/)
    reactant_list(208,0:2) = (/2,18,11/)
    reactant_list(209,0:2) = (/2,18,13/)
    reactant_list(210,0:2) = (/2,18,17/)
    reactant_list(211,0:2) = (/2,18,17/)
    reactant_list(212,0:2) = (/2,35,11/)
    reactant_list(213,0:2) = (/2,35,13/)
    reactant_list(214,0:3) = (/3,2,4,4/)
    reactant_list(215,0:3) = (/3,10,4,4/)
    reactant_list(216,0:3) = (/3,10,4,6/)
    reactant_list(217,0:3) = (/3,12,4,6/)
    reactant_list(218,0:3) = (/3,1,1,4/)

    ! product list ----------------------------
    product_list(1,0:2) = (/2,2,3/)
    product_list(2,0:2) = (/2,5,3/)
    product_list(3,0:3) = (/3,2,3,1/)
    product_list(4,0:2) = (/2,7,3/)
    product_list(5,0:2) = (/2,9,3/)
    product_list(6,0:3) = (/3,10,3,1/)
    product_list(7,0:2) = (/2,5,3/)
    product_list(8,0:2) = (/2,12,3/)
    product_list(9,0:2) = (/2,14,3/)
    product_list(10,0:3) = (/3,15,3,1/)
    product_list(11,0:2) = (/2,12,3/)
    product_list(12,0:2) = (/2,16,3/)
    product_list(13,0:2) = (/2,18,3/)
    product_list(14,0:3) = (/3,19,3,1/)
    product_list(15,0:2) = (/2,14,3/)
    product_list(16,0:2) = (/2,15,3/)
    product_list(17,0:2) = (/2,12,3/)
    product_list(18,0:1) = (/1,1/)
    product_list(19,0:1) = (/1,6/)
    product_list(20,0:2) = (/2,1,6/)
    product_list(21,0:2) = (/2,1,1/)
    product_list(22,0:2) = (/2,4,1/)
    product_list(23,0:3) = (/3,1,1,1/)
    product_list(24,0:1) = (/1,23/)
    product_list(25,0:2) = (/2,23,1/)
    product_list(26,0:2) = (/2,26,1/)
    product_list(27,0:2) = (/2,27,1/)
    product_list(28,0:2) = (/2,28,1/)
    product_list(29,0:3) = (/3,27,1,1/)
    product_list(30,0:3) = (/3,27,1,4/)
    product_list(31,0:3) = (/3,28,1,1/)
    product_list(32,0:2) = (/2,23,23/)
    product_list(33,0:2) = (/2,31,1/)
    product_list(34,0:2) = (/2,26,23/)
    product_list(35,0:2) = (/2,32,1/)
    product_list(36,0:2) = (/2,26,26/)
    product_list(37,0:2) = (/2,11,1/)
    product_list(38,0:2) = (/2,27,26/)
    product_list(39,0:2) = (/2,33,1/)
    product_list(40,0:2) = (/2,27,27/)
    product_list(41,0:2) = (/2,13,1/)
    product_list(42,0:2) = (/2,28,27/)
    product_list(43,0:2) = (/2,34,1/)
    product_list(44,0:2) = (/2,28,28/)
    product_list(45,0:2) = (/2,17,1/)
    product_list(46,0:1) = (/1,0/)
    product_list(47,0:1) = (/1,0/)
    product_list(48,0:1) = (/1,0/)
    product_list(49,0:1) = (/1,0/)
    product_list(50,0:1) = (/1,0/)
    product_list(51,0:1) = (/1,0/)
    product_list(52,0:1) = (/1,0/)
    product_list(53,0:1) = (/1,0/)
    product_list(54,0:1) = (/1,0/)
    product_list(55,0:1) = (/1,0/)
    product_list(56,0:1) = (/1,0/)
    product_list(57,0:1) = (/1,0/)
    product_list(58,0:1) = (/1,0/)
    product_list(59,0:1) = (/1,0/)
    product_list(60,0:1) = (/1,0/)
    product_list(61,0:2) = (/2,21,1/)
    product_list(62,0:2) = (/2,2,4/)
    product_list(63,0:2) = (/2,20,1/)
    product_list(64,0:2) = (/2,29,1/)
    product_list(65,0:2) = (/2,9,4/)
    product_list(66,0:3) = (/3,10,1,4/)
    product_list(67,0:2) = (/2,15,1/)
    product_list(68,0:2) = (/2,12,4/)
    product_list(69,0:2) = (/2,18,4/)
    product_list(70,0:3) = (/3,19,1,4/)
    product_list(71,0:3) = (/3,14,4,4/)
    product_list(72,0:4) = (/4,15,1,4,4/)
    product_list(73,0:4) = (/4,12,4,4,4/)
    product_list(74,0:2) = (/2,5,6/)
    product_list(75,0:3) = (/3,2,1,6/)
    product_list(76,0:3) = (/3,2,1,6/)
    product_list(77,0:3) = (/3,2,28,6/)
    product_list(78,0:4) = (/4,24,4,1,6/)
    product_list(79,0:3) = (/3,25,4,6/)
    product_list(80,0:3) = (/3,10,1,6/)
    product_list(81,0:2) = (/2,9,6/)
    product_list(82,0:2) = (/2,12,6/)
    product_list(83,0:3) = (/3,16,1,6/)
    product_list(84,0:3) = (/3,30,4,6/)
    product_list(85,0:3) = (/3,24,26,6/)
    product_list(86,0:2) = (/2,14,6/)
    product_list(87,0:3) = (/3,15,1,6/)
    product_list(88,0:3) = (/3,12,4,6/)
    product_list(89,0:4) = (/4,16,1,4,6/)
    product_list(90,0:3) = (/3,25,27,6/)
    product_list(91,0:3) = (/3,14,4,6/)
    product_list(92,0:4) = (/4,15,1,4,6/)
    product_list(93,0:4) = (/4,12,4,4,6/)
    product_list(94,0:2) = (/2,10,4/)
    product_list(95,0:2) = (/2,9,1/)
    product_list(96,0:3) = (/3,15,4,4/)
    product_list(97,0:3) = (/3,14,1,4/)
    product_list(98,0:2) = (/2,19,4/)
    product_list(99,0:2) = (/2,5,1/)
    product_list(100,0:2) = (/2,21,6/)
    product_list(101,0:2) = (/2,5,6/)
    product_list(102,0:3) = (/3,14,1,6/)
    product_list(103,0:3) = (/3,15,4,6/)
    product_list(104,0:4) = (/4,15,4,4,6/)
    product_list(105,0:3) = (/3,19,4,6/)
    product_list(106,0:2) = (/2,29,4/)
    product_list(107,0:2) = (/2,15,4/)
    product_list(108,0:2) = (/2,19,4/)
    product_list(109,0:3) = (/3,15,4,4/)
    product_list(110,0:3) = (/3,19,4,4/)
    product_list(111,0:2) = (/2,12,4/)
    product_list(112,0:2) = (/2,15,1/)
    product_list(113,0:2) = (/2,36,1/)
    product_list(114,0:2) = (/2,15,26/)
    product_list(115,0:2) = (/2,14,23/)
    product_list(116,0:3) = (/3,36,1,4/)
    product_list(117,0:2) = (/2,37,4/)
    product_list(118,0:2) = (/2,38,1/)
    product_list(119,0:2) = (/2,12,8/)
    product_list(120,0:2) = (/2,15,28/)
    product_list(121,0:2) = (/2,19,26/)
    product_list(122,0:3) = (/3,38,1,4/)
    product_list(123,0:2) = (/2,22,4/)
    product_list(124,0:2) = (/2,25,1/)
    product_list(125,0:2) = (/2,14,1/)
    product_list(126,0:2) = (/2,15,4/)
    product_list(127,0:3) = (/3,12,1,4/)
    product_list(128,0:2) = (/2,37,1/)
    product_list(129,0:2) = (/2,10,1/)
    product_list(130,0:2) = (/2,19,1/)
    product_list(131,0:2) = (/2,14,4/)
    product_list(132,0:2) = (/2,38,1/)
    product_list(133,0:2) = (/2,19,4/)
    product_list(134,0:2) = (/2,38,4/)
    product_list(135,0:2) = (/2,15,8/)
    product_list(136,0:3) = (/3,38,4,4/)
    product_list(137,0:2) = (/2,40,4/)
    product_list(138,0:2) = (/2,19,8/)
    product_list(139,0:3) = (/3,40,4,4/)
    product_list(140,0:2) = (/2,42,4/)
    product_list(141,0:2) = (/2,29,1/)
    product_list(142,0:2) = (/2,29,28/)
    product_list(143,0:2) = (/2,12,8/)
    product_list(144,0:2) = (/2,15,28/)
    product_list(145,0:3) = (/3,38,1,4/)
    product_list(146,0:2) = (/2,14,8/)
    product_list(147,0:2) = (/2,19,28/)
    product_list(148,0:3) = (/3,40,1,4/)
    product_list(149,0:3) = (/3,14,8,4/)
    product_list(150,0:2) = (/2,9,4/)
    product_list(151,0:2) = (/2,15,8/)
    product_list(152,0:2) = (/2,19,8/)
    product_list(153,0:3) = (/3,19,8,4/)
    product_list(154,0:2) = (/2,35,8/)
    product_list(155,0:2) = (/2,16,1/)
    product_list(156,0:2) = (/2,16,28/)
    product_list(157,0:2) = (/2,12,27/)
    product_list(158,0:3) = (/3,36,1,4/)
    product_list(159,0:2) = (/2,37,4/)
    product_list(160,0:2) = (/2,38,1/)
    product_list(161,0:2) = (/2,45,1/)
    product_list(162,0:1) = (/1,0/)
    product_list(163,0:2) = (/2,12,1/)
    product_list(164,0:2) = (/2,12,28/)
    product_list(165,0:2) = (/2,38,4/)
    product_list(166,0:2) = (/2,39,1/)
    product_list(167,0:1) = (/1,40/)
    product_list(168,0:2) = (/2,46,1/)
    product_list(169,0:1) = (/1,0/)
    product_list(170,0:2) = (/2,15,1/)
    product_list(171,0:2) = (/2,39,4/)
    product_list(172,0:2) = (/2,40,1/)
    product_list(173,0:2) = (/2,46,4/)
    product_list(174,0:2) = (/2,47,1/)
    product_list(175,0:2) = (/2,14,11/)
    product_list(176,0:2) = (/2,38,28/)
    product_list(177,0:2) = (/2,48,1/)
    product_list(178,0:2) = (/2,14,13/)
    product_list(179,0:2) = (/2,19,33/)
    product_list(180,0:3) = (/3,38,28,4/)
    product_list(181,0:2) = (/2,39,8/)
    product_list(182,0:2) = (/2,40,28/)
    product_list(183,0:3) = (/3,48,1,4/)
    product_list(184,0:2) = (/2,49,1/)
    product_list(185,0:2) = (/2,12,4/)
    product_list(186,0:2) = (/2,40,4/)
    product_list(187,0:2) = (/2,47,4/)
    product_list(188,0:2) = (/2,19,11/)
    product_list(189,0:2) = (/2,19,13/)
    product_list(190,0:2) = (/2,40,8/)
    product_list(191,0:2) = (/2,49,4/)
    product_list(192,0:2) = (/2,15,4/)
    product_list(193,0:2) = (/2,38,28/)
    product_list(194,0:2) = (/2,48,1/)
    product_list(195,0:2) = (/2,40,28/)
    product_list(196,0:2) = (/2,49,1/)
    product_list(197,0:2) = (/2,41,8/)
    product_list(198,0:2) = (/2,42,28/)
    product_list(199,0:2) = (/2,14,4/)
    product_list(200,0:2) = (/2,42,4/)
    product_list(201,0:2) = (/2,38,8/)
    product_list(202,0:2) = (/2,48,4/)
    product_list(203,0:2) = (/2,40,8/)
    product_list(204,0:2) = (/2,50,4/)
    product_list(205,0:2) = (/2,19,4/)
    product_list(206,0:2) = (/2,19,33/)
    product_list(207,0:2) = (/2,40,28/)
    product_list(208,0:2) = (/2,49,1/)
    product_list(209,0:2) = (/2,14,17/)
    product_list(210,0:2) = (/2,43,8/)
    product_list(211,0:2) = (/2,44,28/)
    product_list(212,0:2) = (/2,15,17/)
    product_list(213,0:2) = (/2,19,17/)
    product_list(214,0:2) = (/2,21,4/)
    product_list(215,0:2) = (/2,29,4/)
    product_list(216,0:2) = (/2,29,6/)
    product_list(217,0:2) = (/2,14,6/)
    product_list(218,0:2) = (/2,4,4/)

  end subroutine get_reactant_product_list


  subroutine get_reaction_type_list(reaction_type)
    implicit none
    integer, parameter :: nch = 218
    integer, intent(out) :: reaction_type(1:nch)

    ! reaction type
    reaction_type = 0
    reaction_type(1) = 1
    reaction_type(2) = 1
    reaction_type(3) = 1
    reaction_type(4) = 1
    reaction_type(5) = 1
    reaction_type(6) = 1
    reaction_type(7) = 1
    reaction_type(8) = 1
    reaction_type(9) = 1
    reaction_type(10) = 1
    reaction_type(11) = 1
    reaction_type(12) = 1
    reaction_type(13) = 1
    reaction_type(14) = 1
    reaction_type(15) = 1
    reaction_type(16) = 1
    reaction_type(17) = 1

  end subroutine get_reaction_type_list


  integer function get_number_of_wavelength_bins()
    get_number_of_wavelength_bins = 1000
  end function


  subroutine get_wavelength_bin(lambda, dlambda)
    implicit none
    real(dp), intent(inout) :: lambda(1:), dlambda(1:)
    integer i
    do i = 1, 1000
      dlambda(i) = 1.0_dp
      lambda(i) = 0.0_dp + 0.5_dp * dlambda(i) + dble(i-1) * 1.0_dp
    end do
  end subroutine get_wavelength_bin


  subroutine get_solar_flux(lambda, dlambda, solar_flux)
    implicit none
    real(dp), intent(in) :: lambda(1:), dlambda(1:)
    real(dp), intent(out) :: solar_flux(1:)
    integer, parameter :: nwl = 1000
    character(len=256) fname, unit1, unit2

    fname = './UV/ref_solar_irradiance_whi-2008_ver2_1.dat'
    unit1 = 'nm'
    unit2 = 'W/m^2/nm'

    call get_solar_flux_data(lambda, dlambda, nwl, fname, unit1, unit2, & ! in
      &                      solar_flux                                 ) ! out

  end subroutine get_solar_flux


  subroutine p__PROTEUS_source(var_species_list, fix_species_list, & ! in:  name of variable species and fixed species
    &                          nsp_var_out, nsp_fix_out,           & ! in:  number of variable species and fixed species
    &                          nx, ny, nz,                         & ! in:  number of grids in x, y, z direction
    &                          T_n, T_i, T_e,                      & ! in:  temperature of neutrals, ions and electrons
    &                          v_var, v_fix,                       & ! in:  three dimensional velocity vector of variable and fixed species
    &                          n_var, n_fix,                       & ! in:  number density of variable and fixed species
    &                          J_rate_in,                          & ! in:  Reaction rate input, e.g., photolysis rate, rate calculated by other model
    &                          prod, loss, k_coef                  ) ! out: production and loss rates for variable species, and reaction rate coefficients
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    !   * Both arrays are strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables and fixed, respectively. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of species_list string arrays.
    !
    ! - nsp_var_out, nsp_fix_out
    !   * The number of chemical species to be treated as variables and fixed, respectively. 
    !
    ! - nx, ny, nz
    !   * The number of grids in x, y, z direction. 
    !   * Please input "1, 1, nz" if the simulation is vertical 1D, and "1, 1, 1" if it is 0D.
    !
    ! - T_n(nx,ny,nz), T_i(nx,ny,nz,nsp_fix_out+nsp_var_out), T_e(nx,ny,nz)
    !   * T_n, T_i, and T_e are the neutral, ion and electron temperatures in [K] in the simulation grids. 
    !   * Ion temperature Ti has the number of elements equal to the number of variable species and simulation grids,  
    !     and should be given for each variable species ALSO FOR NEUTRALS AND ELECTRONS IF THEY ARE VARIABLES  
    !     because PROTEUS does not give a dedicated arrays only for ions to avoid the complexity of the code interface.
    !   * You can simply define "0.0 [K]" for neutral and electrons in T_i array. 
    !
    ! - v_var(3,nx,ny,nz,nsp_var_out), v_fix(3,nx,ny,nz,nsp_fix_out)
    !   * Three dimensional velocity vector of variable and fixed species in [cm s^-1] in simulation grids. 
    !     v_*(1,ix,iy,iz,i) = x component of the velocity of ith species at (ix,iy,iz) grid. 
    !     v_*(2,ix,iy,iz,i) = y component of the velocity of ith species at (ix,iy,iz) grid. 
    !     v_*(3,ix,iy,iz,i) = z component of the velocity of ith species at (ix,iy,iz) grid. 
    !
    ! - n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    !   * Number density of variable and fixed species in [cm^-3] in simulation grids. 
    !
    ! - J_rate_in(nx,ny,nz,nch)
    !   * Users can input reaction rate coefficient calculated by other modules, such as photolysis rate. Units are arbitrary.
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - prod(nx,ny,nz,nsp_var_out), loss(nx,ny,nz,nsp_var_out)
    !   * Production and loss rates for variable species in [cm^-3 s^-1]. 
    !
    ! - k_coef(nx,ny,nz,nch)
    !   * Reaction rate coefficients for each reaction in the simulation grids in the units as follows:
    !     [s^-1] for reactions with only one reactant 
    !     [cm^3 s^-1] for two-body reactions
    !     [cm^6 s^-1] for three-body reactions
    !
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer, parameter :: nsp = 52
    integer, parameter :: nsp_fix_in = 17
    integer, parameter :: nsp_var_in = 35
    integer, parameter :: nch = 218
    integer,          intent(in)  :: nsp_var_out, nsp_fix_out, nx, ny, nz 
    character(len=*), intent(in)  :: var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    real(dp),         intent(in)  :: T_n(nx,ny,nz), T_i(nx,ny,nz,nsp_fix_out+nsp_var_out), T_e(nx,ny,nz)
    real(dp),         intent(in)  :: v_var(3,nx,ny,nz,nsp_var_out), v_fix(3,nx,ny,nz,nsp_fix_out)
    real(dp),         intent(in)  :: n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    real(dp),         intent(in)  :: J_rate_in(nx,ny,nz,nch)
    real(dp),         intent(out) :: prod(nx,ny,nz,nsp_var_out), loss(nx,ny,nz,nsp_var_out), k_coef(nx,ny,nz,nch)
    integer isp, jsp, ksp 
    integer ix, iy, iz 
    integer  all_in_var_out(nsp), var_out_all_in(nsp_var_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_fix_out(nsp), fix_out_all_in(nsp_fix_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_var_in(nsp) ! 
    real(dp) P(nsp_var_in), L(nsp_var_in), k(nch), k0, kinf, k2, k3, n(nsp), n_tot, Tn, Te, Ti(nsp), T_tmp, v(3,nsp)
    character(len = 256) species(nsp)

    ! species list ----------------------------
    ! 
    ! H, H+, e-, H2, H2+, He, He+, CH4, CH4+, CH3+, C2H2, C2H2+, C2H4, C2H4+, C2H3+, C2H+, C2H6, C2H6+, C2H5+, HeH+, H3+, C+, C, CH+, CH2+, CH, CH2, CH3, CH5+, C2+, C2, C2H, C2H3, C2H5, C2H7+, C3H+, C3H2+, C3H3+, C3H4+, C3H5+, C3H6+, C3H7+, C3H8+, C3H9+, C4H+, C4H2+, C4H3+, C4H5+, C4H7+, C4H9+, H2(v2), H2(v4)
    species(1) = "H"
    species(2) = "H+"
    species(3) = "e-"
    species(4) = "H2"
    species(5) = "H2+"
    species(6) = "He"
    species(7) = "He+"
    species(8) = "CH4"
    species(9) = "CH4+"
    species(10) = "CH3+"
    species(11) = "C2H2"
    species(12) = "C2H2+"
    species(13) = "C2H4"
    species(14) = "C2H4+"
    species(15) = "C2H3+"
    species(16) = "C2H+"
    species(17) = "C2H6"
    species(18) = "C2H6+"
    species(19) = "C2H5+"
    species(20) = "HeH+"
    species(21) = "H3+"
    species(22) = "C+"
    species(23) = "C"
    species(24) = "CH+"
    species(25) = "CH2+"
    species(26) = "CH"
    species(27) = "CH2"
    species(28) = "CH3"
    species(29) = "CH5+"
    species(30) = "C2+"
    species(31) = "C2"
    species(32) = "C2H"
    species(33) = "C2H3"
    species(34) = "C2H5"
    species(35) = "C2H7+"
    species(36) = "C3H+"
    species(37) = "C3H2+"
    species(38) = "C3H3+"
    species(39) = "C3H4+"
    species(40) = "C3H5+"
    species(41) = "C3H6+"
    species(42) = "C3H7+"
    species(43) = "C3H8+"
    species(44) = "C3H9+"
    species(45) = "C4H+"
    species(46) = "C4H2+"
    species(47) = "C4H3+"
    species(48) = "C4H5+"
    species(49) = "C4H7+"
    species(50) = "C4H9+"
    species(51) = "H2(v2)"
    species(52) = "H2(v4)"

    ! Converting index --------------------------
    fix_out_all_in = 0
    all_in_fix_out = 0
    var_out_all_in = 0
    all_in_var_out = 0
    all_in_var_in  = 0

    all_in_var_in(1) = 1 ! H: variable
    all_in_var_in(2) = 2 ! H+: variable
    all_in_var_in(5) = 3 ! H2+: variable
    all_in_var_in(7) = 4 ! He+: variable
    all_in_var_in(9) = 5 ! CH4+: variable
    all_in_var_in(10) = 6 ! CH3+: variable
    all_in_var_in(12) = 7 ! C2H2+: variable
    all_in_var_in(14) = 8 ! C2H4+: variable
    all_in_var_in(15) = 9 ! C2H3+: variable
    all_in_var_in(16) = 10 ! C2H+: variable
    all_in_var_in(18) = 11 ! C2H6+: variable
    all_in_var_in(19) = 12 ! C2H5+: variable
    all_in_var_in(20) = 13 ! HeH+: variable
    all_in_var_in(21) = 14 ! H3+: variable
    all_in_var_in(22) = 15 ! C+: variable
    all_in_var_in(24) = 16 ! CH+: variable
    all_in_var_in(25) = 17 ! CH2+: variable
    all_in_var_in(29) = 18 ! CH5+: variable
    all_in_var_in(30) = 19 ! C2+: variable
    all_in_var_in(35) = 20 ! C2H7+: variable
    all_in_var_in(36) = 21 ! C3H+: variable
    all_in_var_in(37) = 22 ! C3H2+: variable
    all_in_var_in(38) = 23 ! C3H3+: variable
    all_in_var_in(39) = 24 ! C3H4+: variable
    all_in_var_in(40) = 25 ! C3H5+: variable
    all_in_var_in(41) = 26 ! C3H6+: variable
    all_in_var_in(42) = 27 ! C3H7+: variable
    all_in_var_in(43) = 28 ! C3H8+: variable
    all_in_var_in(44) = 29 ! C3H9+: variable
    all_in_var_in(45) = 30 ! C4H+: variable
    all_in_var_in(46) = 31 ! C4H2+: variable
    all_in_var_in(47) = 32 ! C4H3+: variable
    all_in_var_in(48) = 33 ! C4H5+: variable
    all_in_var_in(49) = 34 ! C4H7+: variable
    all_in_var_in(50) = 35 ! C4H9+: variable

    do isp = 1, nsp_fix_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(fix_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          fix_out_all_in(isp) = jsp 
          all_in_fix_out(jsp) = isp 
        end if 
      end do 
    end do 

    do isp = 1, nsp_var_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(var_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          var_out_all_in(isp) = jsp 
          all_in_var_out(jsp) = isp 
        end if 
      end do 
    end do 

    ! Calculating production and loss rates ---

    prod = 0.0_dp
    loss = 0.0_dp

    do iz = 1, nz
    do iy = 1, ny
    do ix = 1, nx

      Tn = T_n(ix,iy,iz)
      Te = T_e(ix,iy,iz)

      n = 0.0_dp
      n_tot = 0.0_dp
      do isp = 1, nsp_fix_out
        jsp = fix_out_all_in(isp)
        if (jsp/=0) then 
          n(jsp) = n_fix(ix,iy,iz,isp) 
          v(1:3,jsp) = v_fix(1:3,ix,iy,iz,isp) 
          Ti(jsp) = T_i(ix,iy,iz,isp) 
          if (trim(ADJUSTL(species(jsp)))/='M') n_tot = n_tot + n(jsp)
        end if 
      end do 
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        if (jsp/=0) then
          n(jsp) = n_var(ix,iy,iz,isp) 
          v(1:3,jsp) = v_var(1:3,ix,iy,iz,isp) 
          Ti(jsp) = T_i(ix,iy,iz,nsp_fix_out+isp)
          if (trim(ADJUSTL(species(jsp)))/='M') n_tot = n_tot + n(jsp)
        end if 
      end do 

      ! Reaction rate coefficient -----------------
      k = 0.0_dp
      k(1) = J_rate_in(ix,iy,iz,1)
      k(2) = J_rate_in(ix,iy,iz,2)
      k(3) = J_rate_in(ix,iy,iz,3)
      k(4) = J_rate_in(ix,iy,iz,4)
      k(5) = J_rate_in(ix,iy,iz,5)
      k(6) = J_rate_in(ix,iy,iz,6)
      k(7) = J_rate_in(ix,iy,iz,7)
      k(8) = J_rate_in(ix,iy,iz,8)
      k(9) = J_rate_in(ix,iy,iz,9)
      k(10) = J_rate_in(ix,iy,iz,10)
      k(11) = J_rate_in(ix,iy,iz,11)
      k(12) = J_rate_in(ix,iy,iz,12)
      k(13) = J_rate_in(ix,iy,iz,13)
      k(14) = J_rate_in(ix,iy,iz,14)
      k(15) = J_rate_in(ix,iy,iz,15)
      k(16) = J_rate_in(ix,iy,iz,16)
      k(17) = J_rate_in(ix,iy,iz,17)
      k(18) = 4.0e-12_dp*(Te/250.0_dp)**(-0.7_dp)
      k(19) = 4.0e-12_dp*(Te/250.0_dp)**(-0.7_dp)
      k(20) = 1.0e-8_dp*(Te/300.0_dp)**(-0.6_dp)
      k(21) = 2.3e-7_dp*(Te/300.0_dp)**(-0.4_dp)
      k(22) = 4.4e-8_dp*(Te/300.0_dp)**(-0.5_dp)
      k(23) = 5.6e-8_dp*(Te/300.0_dp)**(-0.5_dp)
      k(24) = 4.0e-12_dp*(Te/250.0_dp)**(-0.7_dp)
      k(25) = 1.5e-7_dp*(Te/300.0_dp)**(-0.42_dp)
      k(26) = 2.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(27) = 3.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(28) = 3.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(29) = 3.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(30) = 8.8e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(31) = 2.2e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(32) = 3.0e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(33) = 2.7e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(34) = 2.7e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(35) = 2.7e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(36) = 2.7e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(37) = 4.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(38) = 4.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(39) = 3.0e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(40) = 3.0e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(41) = 7.4e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(42) = 7.4e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(43) = 3.0e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(44) = 3.0e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(45) = 3.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(46) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(47) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(48) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(49) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(50) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(51) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(52) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(53) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(54) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(55) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(56) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(57) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(58) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(59) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(60) = 7.5e-7_dp*(Te/300.0_dp)**(-0.5_dp)
      k(61) = 2.00e-9_dp
      k(62) = 6.40e-10_dp
      k(63) = 1.40e-10_dp
      k(64) = 1.10e-10_dp
      k(65) = 1.41e-9_dp
      k(66) = 2.28e-9_dp
      k(67) = 4.77e-10_dp
      k(68) = 4.82e-9_dp
      k(69) = 2.94e-10_dp
      k(70) = 1.37e-9_dp
      k(71) = 2.35e-9_dp
      k(72) = 6.86e-10_dp
      k(73) = 1.96e-10_dp
      k(74) = 9.35e-15_dp
      k(75) = 4.57e-14_dp
      k(76) = 1.0e-9_dp
      k(77) = 4.76e-10_dp
      k(78) = 2.38e-10_dp
      k(79) = 8.50e-10_dp
      k(80) = 8.50e-11_dp
      k(81) = 5.10e-11_dp
      k(82) = 2.45e-10_dp
      k(83) = 8.75e-10_dp
      k(84) = 1.61e-9_dp
      k(85) = 7.70e-10_dp
      k(86) = 2.38e-10_dp
      k(87) = 1.70e-10_dp
      k(88) = 2.18e-9_dp
      k(89) = 4.42e-10_dp
      k(90) = 4.08e-10_dp
      k(91) = 4.20e-10_dp
      k(92) = 1.74e-9_dp
      k(93) = 8.40e-10_dp
      k(94) = 3.69e-9_dp
      k(95) = 0.81e-9_dp
      k(96) = 1.30e-9_dp
      k(97) = 1.30e-9_dp
      k(98) = 1.30e-9_dp
      k(99) = 2.0e-9_dp
      k(100) = 1.50e-9_dp
      k(101) = 9.10e-10_dp
      k(102) = 7.00e-10_dp
      k(103) = 2.10e-9_dp
      k(104) = 1.05e-9_dp
      k(105) = 1.05e-9_dp
      k(106) = 2.40e-9_dp
      k(107) = 2.90e-9_dp
      k(108) = 0.69e-9_dp
      k(109) = 1.61e-9_dp
      k(110) = 2.40e-9_dp
      k(111) = 3.30e-10_dp
      k(112) = 9.80e-10_dp
      k(113) = 2.80e-9_dp
      k(114) = 8.50e-11_dp
      k(115) = 1.70e-10_dp
      k(116) = 8.50e-11_dp
      k(117) = 3.40e-10_dp
      k(118) = 1.00e-9_dp
      k(119) = 1.70e-10_dp
      k(120) = 5.10e-10_dp
      k(121) = 1.70e-10_dp
      k(122) = 8.50e-10_dp
      k(123) = 7.50e-10_dp
      k(124) = 1.20e-9_dp
      k(125) = 6.50e-11_dp
      k(126) = 1.10e-9_dp
      k(127) = 1.40e-10_dp
      k(128) = 2.40e-9_dp
      k(129) = 1.60e-9_dp
      k(130) = 3.60e-10_dp
      k(131) = 8.40e-10_dp
      k(132) = 2.50e-9_dp
      k(133) = 1.20e-9_dp
      k(134) = 1.20e-9_dp
      k(135) = 3.50e-10_dp
      k(136) = 4.60e-11_dp
      k(137) = 5.20e-10_dp
      k(138) = 1.50e-9_dp
      k(139) = 1.60e-10_dp
      k(140) = 1.00e-10_dp
      k(141) = 3.0e-11_dp
      k(142) = 1.50e-9_dp
      k(143) = 1.13e-9_dp
      k(144) = 1.23e-9_dp
      k(145) = 1.51e-10_dp
      k(146) = 1.38e-9_dp
      k(147) = 4.23e-10_dp
      k(148) = 5.52e-11_dp
      k(149) = 1.91e-9_dp
      k(150) = 1.50e-10_dp
      k(151) = 1.56e-9_dp
      k(152) = 1.50e-9_dp
      k(153) = 2.25e-10_dp
      k(154) = 1.28e-9_dp
      k(155) = 1.20e-9_dp
      k(156) = 2.38e-10_dp
      k(157) = 1.82e-10_dp
      k(158) = 1.96e-10_dp
      k(159) = 5.74e-10_dp
      k(160) = 2.10e-10_dp
      k(161) = 1.20e-9_dp
      k(162) = 1.90e-9_dp
      k(163) = 1.70e-9_dp
      k(164) = 3.74e-10_dp
      k(165) = 3.74e-10_dp
      k(166) = 1.32e-10_dp
      k(167) = 2.20e-10_dp
      k(168) = 1.20e-9_dp
      k(169) = 1.71e-9_dp
      k(170) = 1.80e-12_dp
      k(171) = 1.60e-10_dp
      k(172) = 6.40e-10_dp
      k(173) = 5.16e-10_dp
      k(174) = 5.88e-10_dp
      k(175) = 8.96e-10_dp
      k(176) = 2.80e-10_dp
      k(177) = 2.24e-10_dp
      k(178) = 2.63e-10_dp
      k(179) = 1.31e-10_dp
      k(180) = 8.76e-11_dp
      k(181) = 1.46e-11_dp
      k(182) = 7.88e-10_dp
      k(183) = 7.30e-11_dp
      k(184) = 1.31e-10_dp
      k(185) = 1.00e-10_dp
      k(186) = 2.00e-10_dp
      k(187) = 2.16e-10_dp
      k(188) = 9.30e-10_dp
      k(189) = 2.91e-10_dp
      k(190) = 2.48e-10_dp
      k(191) = 8.06e-11_dp
      k(192) = 3.0e-10_dp
      k(193) = 6.73e-10_dp
      k(194) = 2.37e-10_dp
      k(195) = 7.19e-10_dp
      k(196) = 7.11e-11_dp
      k(197) = 3.71e-13_dp
      k(198) = 4.93e-12_dp
      k(199) = 1.0e-11_dp
      k(200) = 4.00e-15_dp
      k(201) = 6.84e-11_dp
      k(202) = 1.22e-10_dp
      k(203) = 3.90e-10_dp
      k(204) = 4.00e-11_dp
      k(205) = 1.0e-10_dp
      k(206) = 2.22e-10_dp
      k(207) = 8.19e-10_dp
      k(208) = 1.29e-10_dp
      k(209) = 1.15e-9_dp
      k(210) = 7.98e-12_dp
      k(211) = 1.10e-11_dp
      k(212) = 1.0e-9_dp
      k(213) = 1.0e-9_dp
      k(214) = 3.2e-29_dp
      k(215) = 3.3e-28_dp
      k(216) = 1.1e-28_dp
      k(217) = 1.2e-27_dp
      k(218) = 2.7e-31_dp*(Tn/1.0_dp)**(-0.6_dp)

      k_coef(ix,iy,iz,:) = k(:)

      ! Production rate ---------------------------
      P = 0.0_dp
      P(1) = k(3)*n(4) &
             & + k(6)*n(8) &
             & + k(10)*n(13) &
             & + k(14)*n(17) &
             & + k(18)*n(2)*n(3) &
             & + k(20)*n(20)*n(3) &
             & + 2.0_dp*k(21)*n(5)*n(3) &
             & + k(22)*n(21)*n(3) &
             & + 3.0_dp*k(23)*n(21)*n(3) &
             & + k(25)*n(24)*n(3) &
             & + k(26)*n(25)*n(3) &
             & + k(27)*n(10)*n(3) &
             & + k(28)*n(9)*n(3) &
             & + 2.0_dp*k(29)*n(9)*n(3) &
             & + k(30)*n(29)*n(3) &
             & + 2.0_dp*k(31)*n(29)*n(3) &
             & + k(33)*n(16)*n(3) &
             & + k(35)*n(12)*n(3) &
             & + k(37)*n(15)*n(3) &
             & + k(39)*n(14)*n(3) &
             & + k(41)*n(19)*n(3) &
             & + k(43)*n(18)*n(3) &
             & + k(45)*n(35)*n(3) &
             & + k(61)*n(5)*n(4) &
             & + k(63)*n(5)*n(6) &
             & + k(64)*n(5)*n(8) &
             & + k(66)*n(5)*n(8) &
             & + k(67)*n(5)*n(11) &
             & + k(70)*n(5)*n(17) &
             & + k(72)*n(5)*n(17) &
             & + k(75)*n(7)*n(4) &
             & + k(76)*n(7)*n(51) &
             & + k(78)*n(7)*n(8) &
             & + k(80)*n(7)*n(8) &
             & + k(83)*n(7)*n(11) &
             & + k(87)*n(7)*n(13) &
             & + k(89)*n(7)*n(13) &
             & + k(92)*n(7)*n(17) &
             & + k(95)*n(2)*n(8) &
             & + k(97)*n(2)*n(17) &
             & + k(99)*n(2)*n(52) &
             & + k(102)*n(20)*n(13) &
             & + k(112)*n(22)*n(8) &
             & + k(113)*n(22)*n(11) &
             & + k(116)*n(22)*n(13) &
             & + k(118)*n(22)*n(13) &
             & + k(122)*n(22)*n(17) &
             & + k(124)*n(24)*n(4) &
             & + k(125)*n(24)*n(8) &
             & + k(127)*n(24)*n(8) &
             & + k(128)*n(24)*n(11) &
             & + k(129)*n(25)*n(4) &
             & + k(130)*n(25)*n(8) &
             & + k(132)*n(25)*n(11) &
             & + k(141)*n(9)*n(4) &
             & + k(145)*n(9)*n(11) &
             & + k(148)*n(9)*n(13) &
             & + k(155)*n(30)*n(4) &
             & + k(158)*n(30)*n(8) &
             & + k(160)*n(30)*n(8) &
             & + k(161)*n(30)*n(11) &
             & + k(163)*n(16)*n(4) &
             & + k(166)*n(16)*n(8) &
             & + k(168)*n(16)*n(11) &
             & + k(170)*n(12)*n(4) &
             & + k(172)*n(12)*n(8) &
             & + k(174)*n(12)*n(11) &
             & + k(177)*n(12)*n(13) &
             & + k(183)*n(12)*n(17) &
             & + k(184)*n(12)*n(17) &
             & + k(194)*n(14)*n(11) &
             & + k(196)*n(14)*n(13) &
             & + k(208)*n(18)*n(11)
      P(2) = k(1)*n(1) &
             & + k(3)*n(4) &
             & + k(62)*n(5)*n(1) &
             & + k(75)*n(7)*n(4) &
             & + k(76)*n(7)*n(51) &
             & + k(77)*n(7)*n(8)
      P(3) = k(2)*n(4) &
             & + k(7)*n(8) &
             & + k(74)*n(7)*n(4) &
             & + k(99)*n(2)*n(52) &
             & + k(101)*n(20)*n(1)
      P(4) = k(4)*n(6)
      P(5) = k(5)*n(8) &
             & + k(65)*n(5)*n(8) &
             & + k(81)*n(7)*n(8) &
             & + k(95)*n(2)*n(8) &
             & + k(150)*n(29)*n(1)
      P(6) = k(6)*n(8) &
             & + k(66)*n(5)*n(8) &
             & + k(80)*n(7)*n(8) &
             & + k(94)*n(2)*n(8) &
             & + k(129)*n(25)*n(4)
      P(7) = k(8)*n(11) &
             & + k(11)*n(13) &
             & + k(17)*n(17) &
             & + k(68)*n(5)*n(11) &
             & + k(73)*n(5)*n(17) &
             & + k(82)*n(7)*n(11) &
             & + k(88)*n(7)*n(13) &
             & + k(93)*n(7)*n(17) &
             & + k(111)*n(22)*n(8) &
             & + k(119)*n(22)*n(17) &
             & + k(127)*n(24)*n(8) &
             & + k(143)*n(9)*n(11) &
             & + k(157)*n(30)*n(8) &
             & + k(163)*n(16)*n(4) &
             & + k(164)*n(16)*n(8) &
             & + k(185)*n(15)*n(1)
      P(8) = k(9)*n(13) &
             & + k(15)*n(17) &
             & + k(71)*n(5)*n(17) &
             & + k(86)*n(7)*n(13) &
             & + k(91)*n(7)*n(17) &
             & + k(97)*n(2)*n(17) &
             & + k(102)*n(20)*n(13) &
             & + k(115)*n(22)*n(13) &
             & + k(125)*n(24)*n(8) &
             & + k(131)*n(25)*n(8) &
             & + k(146)*n(9)*n(13) &
             & + k(149)*n(9)*n(17) &
             & + k(175)*n(12)*n(13) &
             & + k(178)*n(12)*n(17) &
             & + k(199)*n(19)*n(1) &
             & + k(209)*n(18)*n(13) &
             & + k(217)*n(12)*n(4)*n(6)
      P(9) = k(10)*n(13) &
             & + k(16)*n(17) &
             & + k(67)*n(5)*n(11) &
             & + k(72)*n(5)*n(17) &
             & + k(87)*n(7)*n(13) &
             & + k(92)*n(7)*n(17) &
             & + k(96)*n(2)*n(17) &
             & + k(103)*n(20)*n(13) &
             & + k(104)*n(20)*n(17) &
             & + k(107)*n(21)*n(11) &
             & + k(109)*n(21)*n(13) &
             & + k(112)*n(22)*n(8) &
             & + k(114)*n(22)*n(13) &
             & + k(120)*n(22)*n(17) &
             & + k(126)*n(24)*n(8) &
             & + k(135)*n(10)*n(13) &
             & + k(144)*n(9)*n(11) &
             & + k(151)*n(29)*n(11) &
             & + k(170)*n(12)*n(4) &
             & + k(192)*n(14)*n(1) &
             & + k(212)*n(35)*n(11)
      P(10) = k(12)*n(13) &
             & + k(83)*n(7)*n(11) &
             & + k(89)*n(7)*n(13) &
             & + k(155)*n(30)*n(4) &
             & + k(156)*n(30)*n(8)
      P(11) = k(13)*n(17) &
             & + k(69)*n(5)*n(17)
      P(12) = k(14)*n(17) &
             & + k(70)*n(5)*n(17) &
             & + k(98)*n(2)*n(17) &
             & + k(105)*n(20)*n(17) &
             & + k(108)*n(21)*n(13) &
             & + k(110)*n(21)*n(17) &
             & + k(121)*n(22)*n(17) &
             & + k(130)*n(25)*n(8) &
             & + k(133)*n(10)*n(8) &
             & + k(138)*n(10)*n(17) &
             & + k(147)*n(9)*n(13) &
             & + k(152)*n(29)*n(13) &
             & + k(153)*n(29)*n(17) &
             & + k(179)*n(12)*n(17) &
             & + k(188)*n(15)*n(13) &
             & + k(189)*n(15)*n(17) &
             & + k(205)*n(18)*n(1) &
             & + k(206)*n(18)*n(11) &
             & + k(213)*n(35)*n(13)
      P(13) = k(63)*n(5)*n(6)
      P(14) = k(61)*n(5)*n(4) &
             & + k(100)*n(20)*n(4) &
             & + k(214)*n(2)*n(4)*n(4)
      P(15) = k(123)*n(24)*n(1)
      P(16) = k(78)*n(7)*n(8) &
             & + k(85)*n(7)*n(11)
      P(17) = k(79)*n(7)*n(8) &
             & + k(90)*n(7)*n(13) &
             & + k(124)*n(24)*n(4)
      P(18) = k(64)*n(5)*n(8) &
             & + k(106)*n(21)*n(8) &
             & + k(141)*n(9)*n(4) &
             & + k(142)*n(9)*n(8) &
             & + k(215)*n(10)*n(4)*n(4) &
             & + k(216)*n(10)*n(4)*n(6)
      P(19) = k(84)*n(7)*n(11)
      P(20) = k(154)*n(29)*n(17)
      P(21) = k(113)*n(22)*n(11) &
             & + k(116)*n(22)*n(13) &
             & + k(158)*n(30)*n(8)
      P(22) = k(117)*n(22)*n(13) &
             & + k(128)*n(24)*n(11) &
             & + k(159)*n(30)*n(8)
      P(23) = k(118)*n(22)*n(13) &
             & + k(122)*n(22)*n(17) &
             & + k(132)*n(25)*n(11) &
             & + k(134)*n(10)*n(11) &
             & + k(136)*n(10)*n(13) &
             & + k(145)*n(9)*n(11) &
             & + k(160)*n(30)*n(8) &
             & + k(165)*n(16)*n(8) &
             & + k(176)*n(12)*n(13) &
             & + k(180)*n(12)*n(17) &
             & + k(193)*n(14)*n(11) &
             & + k(201)*n(19)*n(11)
      P(24) = k(166)*n(16)*n(8) &
             & + k(171)*n(12)*n(8) &
             & + k(181)*n(12)*n(17)
      P(25) = k(137)*n(10)*n(13) &
             & + k(139)*n(10)*n(17) &
             & + k(148)*n(9)*n(13) &
             & + k(167)*n(16)*n(8) &
             & + k(172)*n(12)*n(8) &
             & + k(182)*n(12)*n(17) &
             & + k(186)*n(15)*n(8) &
             & + k(190)*n(15)*n(17) &
             & + k(195)*n(14)*n(13) &
             & + k(203)*n(19)*n(13) &
             & + k(207)*n(18)*n(11)
      P(26) = k(197)*n(14)*n(17)
      P(27) = k(140)*n(10)*n(17) &
             & + k(198)*n(14)*n(17) &
             & + k(200)*n(19)*n(8)
      P(28) = k(210)*n(18)*n(17)
      P(29) = k(211)*n(18)*n(17)
      P(30) = k(161)*n(30)*n(11)
      P(31) = k(168)*n(16)*n(11) &
             & + k(173)*n(12)*n(11)
      P(32) = k(174)*n(12)*n(11) &
             & + k(187)*n(15)*n(11)
      P(33) = k(177)*n(12)*n(13) &
             & + k(183)*n(12)*n(17) &
             & + k(194)*n(14)*n(11) &
             & + k(202)*n(19)*n(11)
      P(34) = k(184)*n(12)*n(17) &
             & + k(191)*n(15)*n(17) &
             & + k(196)*n(14)*n(13) &
             & + k(208)*n(18)*n(11)
      P(35) = k(204)*n(19)*n(17)

      ! Loss rate ---------------------------
      L = 0.0_dp
      L(1) = k(1)*n(1) &
             & + k(62)*n(5)*n(1) &
             & + k(101)*n(20)*n(1) &
             & + k(123)*n(24)*n(1) &
             & + k(150)*n(29)*n(1) &
             & + k(185)*n(15)*n(1) &
             & + k(192)*n(14)*n(1) &
             & + k(199)*n(19)*n(1) &
             & + k(205)*n(18)*n(1) &
             & + 2.0_dp*k(218)*n(1)*n(1)*n(4)
      L(2) = k(18)*n(2)*n(3) &
             & + k(94)*n(2)*n(8) &
             & + k(95)*n(2)*n(8) &
             & + k(96)*n(2)*n(17) &
             & + k(97)*n(2)*n(17) &
             & + k(98)*n(2)*n(17) &
             & + k(99)*n(2)*n(52) &
             & + k(214)*n(2)*n(4)*n(4)
      L(3) = k(21)*n(5)*n(3) &
             & + k(61)*n(5)*n(4) &
             & + k(62)*n(5)*n(1) &
             & + k(63)*n(5)*n(6) &
             & + k(64)*n(5)*n(8) &
             & + k(65)*n(5)*n(8) &
             & + k(66)*n(5)*n(8) &
             & + k(67)*n(5)*n(11) &
             & + k(68)*n(5)*n(11) &
             & + k(69)*n(5)*n(17) &
             & + k(70)*n(5)*n(17) &
             & + k(71)*n(5)*n(17) &
             & + k(72)*n(5)*n(17) &
             & + k(73)*n(5)*n(17)
      L(4) = k(19)*n(7)*n(3) &
             & + k(74)*n(7)*n(4) &
             & + k(75)*n(7)*n(4) &
             & + k(76)*n(7)*n(51) &
             & + k(77)*n(7)*n(8) &
             & + k(78)*n(7)*n(8) &
             & + k(79)*n(7)*n(8) &
             & + k(80)*n(7)*n(8) &
             & + k(81)*n(7)*n(8) &
             & + k(82)*n(7)*n(11) &
             & + k(83)*n(7)*n(11) &
             & + k(84)*n(7)*n(11) &
             & + k(85)*n(7)*n(11) &
             & + k(86)*n(7)*n(13) &
             & + k(87)*n(7)*n(13) &
             & + k(88)*n(7)*n(13) &
             & + k(89)*n(7)*n(13) &
             & + k(90)*n(7)*n(13) &
             & + k(91)*n(7)*n(17) &
             & + k(92)*n(7)*n(17) &
             & + k(93)*n(7)*n(17)
      L(5) = k(28)*n(9)*n(3) &
             & + k(29)*n(9)*n(3) &
             & + k(141)*n(9)*n(4) &
             & + k(142)*n(9)*n(8) &
             & + k(143)*n(9)*n(11) &
             & + k(144)*n(9)*n(11) &
             & + k(145)*n(9)*n(11) &
             & + k(146)*n(9)*n(13) &
             & + k(147)*n(9)*n(13) &
             & + k(148)*n(9)*n(13) &
             & + k(149)*n(9)*n(17)
      L(6) = k(27)*n(10)*n(3) &
             & + k(133)*n(10)*n(8) &
             & + k(134)*n(10)*n(11) &
             & + k(135)*n(10)*n(13) &
             & + k(136)*n(10)*n(13) &
             & + k(137)*n(10)*n(13) &
             & + k(138)*n(10)*n(17) &
             & + k(139)*n(10)*n(17) &
             & + k(140)*n(10)*n(17) &
             & + k(215)*n(10)*n(4)*n(4) &
             & + k(216)*n(10)*n(4)*n(6)
      L(7) = k(35)*n(12)*n(3) &
             & + k(36)*n(12)*n(3) &
             & + k(170)*n(12)*n(4) &
             & + k(171)*n(12)*n(8) &
             & + k(172)*n(12)*n(8) &
             & + k(173)*n(12)*n(11) &
             & + k(174)*n(12)*n(11) &
             & + k(175)*n(12)*n(13) &
             & + k(176)*n(12)*n(13) &
             & + k(177)*n(12)*n(13) &
             & + k(178)*n(12)*n(17) &
             & + k(179)*n(12)*n(17) &
             & + k(180)*n(12)*n(17) &
             & + k(181)*n(12)*n(17) &
             & + k(182)*n(12)*n(17) &
             & + k(183)*n(12)*n(17) &
             & + k(184)*n(12)*n(17) &
             & + k(217)*n(12)*n(4)*n(6)
      L(8) = k(39)*n(14)*n(3) &
             & + k(40)*n(14)*n(3) &
             & + k(192)*n(14)*n(1) &
             & + k(193)*n(14)*n(11) &
             & + k(194)*n(14)*n(11) &
             & + k(195)*n(14)*n(13) &
             & + k(196)*n(14)*n(13) &
             & + k(197)*n(14)*n(17) &
             & + k(198)*n(14)*n(17)
      L(9) = k(37)*n(15)*n(3) &
             & + k(38)*n(15)*n(3) &
             & + k(185)*n(15)*n(1) &
             & + k(186)*n(15)*n(8) &
             & + k(187)*n(15)*n(11) &
             & + k(188)*n(15)*n(13) &
             & + k(189)*n(15)*n(17) &
             & + k(190)*n(15)*n(17) &
             & + k(191)*n(15)*n(17)
      L(10) = k(33)*n(16)*n(3) &
             & + k(34)*n(16)*n(3) &
             & + k(163)*n(16)*n(4) &
             & + k(164)*n(16)*n(8) &
             & + k(165)*n(16)*n(8) &
             & + k(166)*n(16)*n(8) &
             & + k(167)*n(16)*n(8) &
             & + k(168)*n(16)*n(11) &
             & + k(169)*n(16)*n(13)
      L(11) = k(43)*n(18)*n(3) &
             & + k(44)*n(18)*n(3) &
             & + k(205)*n(18)*n(1) &
             & + k(206)*n(18)*n(11) &
             & + k(207)*n(18)*n(11) &
             & + k(208)*n(18)*n(11) &
             & + k(209)*n(18)*n(13) &
             & + k(210)*n(18)*n(17) &
             & + k(211)*n(18)*n(17)
      L(12) = k(41)*n(19)*n(3) &
             & + k(42)*n(19)*n(3) &
             & + k(199)*n(19)*n(1) &
             & + k(200)*n(19)*n(8) &
             & + k(201)*n(19)*n(11) &
             & + k(202)*n(19)*n(11) &
             & + k(203)*n(19)*n(13) &
             & + k(204)*n(19)*n(17)
      L(13) = k(20)*n(20)*n(3) &
             & + k(100)*n(20)*n(4) &
             & + k(101)*n(20)*n(1) &
             & + k(102)*n(20)*n(13) &
             & + k(103)*n(20)*n(13) &
             & + k(104)*n(20)*n(17) &
             & + k(105)*n(20)*n(17)
      L(14) = k(22)*n(21)*n(3) &
             & + k(23)*n(21)*n(3) &
             & + k(106)*n(21)*n(8) &
             & + k(107)*n(21)*n(11) &
             & + k(108)*n(21)*n(13) &
             & + k(109)*n(21)*n(13) &
             & + k(110)*n(21)*n(17)
      L(15) = k(24)*n(22)*n(3) &
             & + k(111)*n(22)*n(8) &
             & + k(112)*n(22)*n(8) &
             & + k(113)*n(22)*n(11) &
             & + k(114)*n(22)*n(13) &
             & + k(115)*n(22)*n(13) &
             & + k(116)*n(22)*n(13) &
             & + k(117)*n(22)*n(13) &
             & + k(118)*n(22)*n(13) &
             & + k(119)*n(22)*n(17) &
             & + k(120)*n(22)*n(17) &
             & + k(121)*n(22)*n(17) &
             & + k(122)*n(22)*n(17)
      L(16) = k(25)*n(24)*n(3) &
             & + k(123)*n(24)*n(1) &
             & + k(124)*n(24)*n(4) &
             & + k(125)*n(24)*n(8) &
             & + k(126)*n(24)*n(8) &
             & + k(127)*n(24)*n(8) &
             & + k(128)*n(24)*n(11)
      L(17) = k(26)*n(25)*n(3) &
             & + k(129)*n(25)*n(4) &
             & + k(130)*n(25)*n(8) &
             & + k(131)*n(25)*n(8) &
             & + k(132)*n(25)*n(11)
      L(18) = k(30)*n(29)*n(3) &
             & + k(31)*n(29)*n(3) &
             & + k(150)*n(29)*n(1) &
             & + k(151)*n(29)*n(11) &
             & + k(152)*n(29)*n(13) &
             & + k(153)*n(29)*n(17) &
             & + k(154)*n(29)*n(17)
      L(19) = k(32)*n(30)*n(3) &
             & + k(155)*n(30)*n(4) &
             & + k(156)*n(30)*n(8) &
             & + k(157)*n(30)*n(8) &
             & + k(158)*n(30)*n(8) &
             & + k(159)*n(30)*n(8) &
             & + k(160)*n(30)*n(8) &
             & + k(161)*n(30)*n(11) &
             & + k(162)*n(30)*n(13)
      L(20) = k(45)*n(35)*n(3) &
             & + k(212)*n(35)*n(11) &
             & + k(213)*n(35)*n(13)
      L(21) = k(46)*n(36)*n(3)
      L(22) = k(47)*n(37)*n(3)
      L(23) = k(48)*n(38)*n(3)
      L(24) = k(49)*n(39)*n(3)
      L(25) = k(50)*n(40)*n(3)
      L(26) = k(51)*n(41)*n(3)
      L(27) = k(52)*n(42)*n(3)
      L(28) = k(53)*n(43)*n(3)
      L(29) = k(54)*n(44)*n(3)
      L(30) = k(55)*n(45)*n(3)
      L(31) = k(56)*n(46)*n(3)
      L(32) = k(57)*n(47)*n(3)
      L(33) = k(58)*n(48)*n(3)
      L(34) = k(59)*n(49)*n(3)
      L(35) = k(60)*n(50)*n(3)


      ! Converting index --------------------------
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        ksp = all_in_var_in(jsp)
        prod(ix,iy,iz,isp) = P(ksp)
        loss(ix,iy,iz,isp) = L(ksp)
      end do 

    end do 
    end do 
    end do 

  end subroutine p__PROTEUS_source



  subroutine p__PROTEUS_Jacobian(var_species_list, fix_species_list, & ! in:  name of variable species and fixed species
    &                            nsp_var_out, nsp_fix_out,           & ! in:  number of variable species and fixed species
    &                            nz,                                 & ! in:  number of vertical grids
    &                            n_var, n_fix,                       & ! in:  number density of variable and fixed species
    &                            k_coef,                             & ! in:  reaction rate coefficients
    &                            Jmtx                                ) ! out: chemical Jacobian matrix
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    !   * Both arrays are strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables and fixed, respectively. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of species_list string arrays.
    !
    ! - nsp_var_out, nsp_fix_out
    !   * The number of chemical species to be treated as variables and fixed, respectively. 
    !
    ! - nz
    !   * The number of vertical grids. 
    !
    ! - n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    !   * Number density of variable and fixed species in [cm^-3] in simulation grids. 
    !
    ! - k_coef(nx,ny,nz,nch)
    !   * Reaction rate coefficients for each reaction in the simulation grids in the units as follows:
    !     [s^-1] for reactions with only one reactant 
    !     [cm^3 s^-1] for two-body reactions
    !     [cm^6 s^-1] for three-body reactions
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - Jmtx(1:,1:)
    !   * If nz >= 2, Jmtx is the transposed chemical Jacobian matrix with a dimension (2*nsp_var_out+1,nsp_var_out*nz).
    !   * If nz == 1, Jmtx is the chemical Jacobian matrix with a dimension (nsp_var_out,nsp_var_out).
    !
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer, parameter :: nsp = 52
    integer, parameter :: nsp_fix_in = 17
    integer, parameter :: nsp_var_in = 35
    integer, parameter :: nch = 218
    integer,          intent(in)  :: nsp_var_out, nsp_fix_out, nz
    character(len=*), intent(in)  :: var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    real(dp),         intent(in)  :: n_var(nz,nsp_var_out), n_fix(nz,nsp_fix_out)
    real(dp),         intent(in)  :: k_coef(nz,nch)
    real(dp),         intent(out) :: Jmtx(1:,1:)
    integer isp, jsp, i0, i1, j0, j1, i, j, iz
    
    integer  all_in_var_out(nsp), var_out_all_in(nsp_var_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_fix_out(nsp), fix_out_all_in(nsp_fix_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_var_in(nsp) ! 
    real(dp) k(nch), n(nsp), Jmtx_tmp(nsp_var_in,nsp_var_in)
    character(len = 256) species(nsp)
    
    ! species list ----------------------------
    ! 
    ! H, H+, e-, H2, H2+, He, He+, CH4, CH4+, CH3+, C2H2, C2H2+, C2H4, C2H4+, C2H3+, C2H+, C2H6, C2H6+, C2H5+, HeH+, H3+, C+, C, CH+, CH2+, CH, CH2, CH3, CH5+, C2+, C2, C2H, C2H3, C2H5, C2H7+, C3H+, C3H2+, C3H3+, C3H4+, C3H5+, C3H6+, C3H7+, C3H8+, C3H9+, C4H+, C4H2+, C4H3+, C4H5+, C4H7+, C4H9+, H2(v2), H2(v4)
    species(1) = "H"
    species(2) = "H+"
    species(3) = "e-"
    species(4) = "H2"
    species(5) = "H2+"
    species(6) = "He"
    species(7) = "He+"
    species(8) = "CH4"
    species(9) = "CH4+"
    species(10) = "CH3+"
    species(11) = "C2H2"
    species(12) = "C2H2+"
    species(13) = "C2H4"
    species(14) = "C2H4+"
    species(15) = "C2H3+"
    species(16) = "C2H+"
    species(17) = "C2H6"
    species(18) = "C2H6+"
    species(19) = "C2H5+"
    species(20) = "HeH+"
    species(21) = "H3+"
    species(22) = "C+"
    species(23) = "C"
    species(24) = "CH+"
    species(25) = "CH2+"
    species(26) = "CH"
    species(27) = "CH2"
    species(28) = "CH3"
    species(29) = "CH5+"
    species(30) = "C2+"
    species(31) = "C2"
    species(32) = "C2H"
    species(33) = "C2H3"
    species(34) = "C2H5"
    species(35) = "C2H7+"
    species(36) = "C3H+"
    species(37) = "C3H2+"
    species(38) = "C3H3+"
    species(39) = "C3H4+"
    species(40) = "C3H5+"
    species(41) = "C3H6+"
    species(42) = "C3H7+"
    species(43) = "C3H8+"
    species(44) = "C3H9+"
    species(45) = "C4H+"
    species(46) = "C4H2+"
    species(47) = "C4H3+"
    species(48) = "C4H5+"
    species(49) = "C4H7+"
    species(50) = "C4H9+"
    species(51) = "H2(v2)"
    species(52) = "H2(v4)"

    ! Converting index --------------------------
    fix_out_all_in = 0
    all_in_fix_out = 0
    var_out_all_in = 0
    all_in_var_out = 0
    all_in_var_in  = 0

    all_in_var_in(1) = 1 ! H: variable
    all_in_var_in(2) = 2 ! H+: variable
    all_in_var_in(5) = 3 ! H2+: variable
    all_in_var_in(7) = 4 ! He+: variable
    all_in_var_in(9) = 5 ! CH4+: variable
    all_in_var_in(10) = 6 ! CH3+: variable
    all_in_var_in(12) = 7 ! C2H2+: variable
    all_in_var_in(14) = 8 ! C2H4+: variable
    all_in_var_in(15) = 9 ! C2H3+: variable
    all_in_var_in(16) = 10 ! C2H+: variable
    all_in_var_in(18) = 11 ! C2H6+: variable
    all_in_var_in(19) = 12 ! C2H5+: variable
    all_in_var_in(20) = 13 ! HeH+: variable
    all_in_var_in(21) = 14 ! H3+: variable
    all_in_var_in(22) = 15 ! C+: variable
    all_in_var_in(24) = 16 ! CH+: variable
    all_in_var_in(25) = 17 ! CH2+: variable
    all_in_var_in(29) = 18 ! CH5+: variable
    all_in_var_in(30) = 19 ! C2+: variable
    all_in_var_in(35) = 20 ! C2H7+: variable
    all_in_var_in(36) = 21 ! C3H+: variable
    all_in_var_in(37) = 22 ! C3H2+: variable
    all_in_var_in(38) = 23 ! C3H3+: variable
    all_in_var_in(39) = 24 ! C3H4+: variable
    all_in_var_in(40) = 25 ! C3H5+: variable
    all_in_var_in(41) = 26 ! C3H6+: variable
    all_in_var_in(42) = 27 ! C3H7+: variable
    all_in_var_in(43) = 28 ! C3H8+: variable
    all_in_var_in(44) = 29 ! C3H9+: variable
    all_in_var_in(45) = 30 ! C4H+: variable
    all_in_var_in(46) = 31 ! C4H2+: variable
    all_in_var_in(47) = 32 ! C4H3+: variable
    all_in_var_in(48) = 33 ! C4H5+: variable
    all_in_var_in(49) = 34 ! C4H7+: variable
    all_in_var_in(50) = 35 ! C4H9+: variable

    do isp = 1, nsp_fix_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(fix_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          fix_out_all_in(isp) = jsp 
          all_in_fix_out(jsp) = isp 
        end if 
      end do 
    end do 

    do isp = 1, nsp_var_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(var_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          var_out_all_in(isp) = jsp 
          all_in_var_out(jsp) = isp 
        end if 
      end do 
    end do 

    ! Calculating chemical Jacobian matrix ---
    Jmtx = 0.0_dp

    do iz = 1, nz

      n = 0.0_dp
      do isp = 1, nsp_fix_out
        jsp = fix_out_all_in(isp)
        if (jsp/=0) then 
          n(jsp) = n_fix(iz,isp) 
        end if 
      end do 
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        if (jsp/=0) then
          n(jsp) = n_var(iz,isp) 
        end if 
      end do 

      k(1:nch) = k_coef(iz,1:nch)

      Jmtx_tmp = 0.0_dp

      Jmtx_tmp(1,1) = - k(1) &
                  & - k(62)*n(5) &
                  & - k(101)*n(20) &
                  & - k(123)*n(24) &
                  & - k(150)*n(29) &
                  & - k(185)*n(15) &
                  & - k(192)*n(14) &
                  & - k(199)*n(19) &
                  & - k(205)*n(18) &
                  & - 2.0_dp*k(218)*n(4)*2.0_dp*n(1)

      Jmtx_tmp(1,2) = k(18)*n(3) &
                  & + k(95)*n(8) &
                  & + k(97)*n(17) &
                  & + k(99)*n(52)

      Jmtx_tmp(1,3) = 2.0_dp*k(21)*n(3) &
                  & + k(61)*n(4) &
                  & + k(63)*n(6) &
                  & + k(64)*n(8) &
                  & + k(66)*n(8) &
                  & + k(67)*n(11) &
                  & + k(70)*n(17) &
                  & + k(72)*n(17) &
                  & - k(62)*n(1)

      Jmtx_tmp(1,4) = k(75)*n(4) &
                  & + k(76)*n(51) &
                  & + k(78)*n(8) &
                  & + k(80)*n(8) &
                  & + k(83)*n(11) &
                  & + k(87)*n(13) &
                  & + k(89)*n(13) &
                  & + k(92)*n(17)

      Jmtx_tmp(1,5) = k(28)*n(3) &
                  & + 2.0_dp*k(29)*n(3) &
                  & + k(141)*n(4) &
                  & + k(145)*n(11) &
                  & + k(148)*n(13)

      Jmtx_tmp(1,6) = k(27)*n(3)

      Jmtx_tmp(1,7) = k(35)*n(3) &
                  & + k(170)*n(4) &
                  & + k(172)*n(8) &
                  & + k(174)*n(11) &
                  & + k(177)*n(13) &
                  & + k(183)*n(17) &
                  & + k(184)*n(17)

      Jmtx_tmp(1,8) = k(39)*n(3) &
                  & + k(194)*n(11) &
                  & + k(196)*n(13) &
                  & - k(192)*n(1)

      Jmtx_tmp(1,9) = k(37)*n(3) &
                  & - k(185)*n(1)

      Jmtx_tmp(1,10) = k(33)*n(3) &
                  & + k(163)*n(4) &
                  & + k(166)*n(8) &
                  & + k(168)*n(11)

      Jmtx_tmp(1,11) = k(43)*n(3) &
                  & + k(208)*n(11) &
                  & - k(205)*n(1)

      Jmtx_tmp(1,12) = k(41)*n(3) &
                  & - k(199)*n(1)

      Jmtx_tmp(1,13) = k(20)*n(3) &
                  & + k(102)*n(13) &
                  & - k(101)*n(1)

      Jmtx_tmp(1,14) = k(22)*n(3) &
                  & + 3.0_dp*k(23)*n(3)

      Jmtx_tmp(1,15) = k(112)*n(8) &
                  & + k(113)*n(11) &
                  & + k(116)*n(13) &
                  & + k(118)*n(13) &
                  & + k(122)*n(17)

      Jmtx_tmp(1,16) = k(25)*n(3) &
                  & + k(124)*n(4) &
                  & + k(125)*n(8) &
                  & + k(127)*n(8) &
                  & + k(128)*n(11) &
                  & - k(123)*n(1)

      Jmtx_tmp(1,17) = k(26)*n(3) &
                  & + k(129)*n(4) &
                  & + k(130)*n(8) &
                  & + k(132)*n(11)

      Jmtx_tmp(1,18) = k(30)*n(3) &
                  & + 2.0_dp*k(31)*n(3) &
                  & - k(150)*n(1)

      Jmtx_tmp(1,19) = k(155)*n(4) &
                  & + k(158)*n(8) &
                  & + k(160)*n(8) &
                  & + k(161)*n(11)

      Jmtx_tmp(1,20) = k(45)*n(3)

      Jmtx_tmp(2,1) = k(1) &
                  & + k(62)*n(5)

      Jmtx_tmp(2,2) = - k(18)*n(3) &
                  & - k(94)*n(8) &
                  & - k(95)*n(8) &
                  & - k(96)*n(17) &
                  & - k(97)*n(17) &
                  & - k(98)*n(17) &
                  & - k(99)*n(52) &
                  & - k(214)*n(4)*n(4)

      Jmtx_tmp(2,3) = k(62)*n(1)

      Jmtx_tmp(2,4) = k(75)*n(4) &
                  & + k(76)*n(51) &
                  & + k(77)*n(8)

      Jmtx_tmp(3,1) = k(101)*n(20) &
                  & - k(62)*n(5)

      Jmtx_tmp(3,2) = k(99)*n(52)

      Jmtx_tmp(3,3) = - k(21)*n(3) &
                  & - k(61)*n(4) &
                  & - k(62)*n(1) &
                  & - k(63)*n(6) &
                  & - k(64)*n(8) &
                  & - k(65)*n(8) &
                  & - k(66)*n(8) &
                  & - k(67)*n(11) &
                  & - k(68)*n(11) &
                  & - k(69)*n(17) &
                  & - k(70)*n(17) &
                  & - k(71)*n(17) &
                  & - k(72)*n(17) &
                  & - k(73)*n(17)

      Jmtx_tmp(3,4) = k(74)*n(4)

      Jmtx_tmp(3,13) = k(101)*n(1)

      Jmtx_tmp(4,4) = - k(19)*n(3) &
                  & - k(74)*n(4) &
                  & - k(75)*n(4) &
                  & - k(76)*n(51) &
                  & - k(77)*n(8) &
                  & - k(78)*n(8) &
                  & - k(79)*n(8) &
                  & - k(80)*n(8) &
                  & - k(81)*n(8) &
                  & - k(82)*n(11) &
                  & - k(83)*n(11) &
                  & - k(84)*n(11) &
                  & - k(85)*n(11) &
                  & - k(86)*n(13) &
                  & - k(87)*n(13) &
                  & - k(88)*n(13) &
                  & - k(89)*n(13) &
                  & - k(90)*n(13) &
                  & - k(91)*n(17) &
                  & - k(92)*n(17) &
                  & - k(93)*n(17)

      Jmtx_tmp(5,1) = k(150)*n(29)

      Jmtx_tmp(5,2) = k(95)*n(8)

      Jmtx_tmp(5,3) = k(65)*n(8)

      Jmtx_tmp(5,4) = k(81)*n(8)

      Jmtx_tmp(5,5) = - k(28)*n(3) &
                  & - k(29)*n(3) &
                  & - k(141)*n(4) &
                  & - k(142)*n(8) &
                  & - k(143)*n(11) &
                  & - k(144)*n(11) &
                  & - k(145)*n(11) &
                  & - k(146)*n(13) &
                  & - k(147)*n(13) &
                  & - k(148)*n(13) &
                  & - k(149)*n(17)

      Jmtx_tmp(5,18) = k(150)*n(1)

      Jmtx_tmp(6,2) = k(94)*n(8)

      Jmtx_tmp(6,3) = k(66)*n(8)

      Jmtx_tmp(6,4) = k(80)*n(8)

      Jmtx_tmp(6,6) = - k(27)*n(3) &
                  & - k(133)*n(8) &
                  & - k(134)*n(11) &
                  & - k(135)*n(13) &
                  & - k(136)*n(13) &
                  & - k(137)*n(13) &
                  & - k(138)*n(17) &
                  & - k(139)*n(17) &
                  & - k(140)*n(17) &
                  & - k(215)*n(4)*n(4) &
                  & - k(216)*n(4)*n(6)

      Jmtx_tmp(6,17) = k(129)*n(4)

      Jmtx_tmp(7,1) = k(185)*n(15)

      Jmtx_tmp(7,3) = k(68)*n(11) &
                  & + k(73)*n(17)

      Jmtx_tmp(7,4) = k(82)*n(11) &
                  & + k(88)*n(13) &
                  & + k(93)*n(17)

      Jmtx_tmp(7,5) = k(143)*n(11)

      Jmtx_tmp(7,7) = - k(35)*n(3) &
                  & - k(36)*n(3) &
                  & - k(170)*n(4) &
                  & - k(171)*n(8) &
                  & - k(172)*n(8) &
                  & - k(173)*n(11) &
                  & - k(174)*n(11) &
                  & - k(175)*n(13) &
                  & - k(176)*n(13) &
                  & - k(177)*n(13) &
                  & - k(178)*n(17) &
                  & - k(179)*n(17) &
                  & - k(180)*n(17) &
                  & - k(181)*n(17) &
                  & - k(182)*n(17) &
                  & - k(183)*n(17) &
                  & - k(184)*n(17) &
                  & - k(217)*n(4)*n(6)

      Jmtx_tmp(7,9) = k(185)*n(1)

      Jmtx_tmp(7,10) = k(163)*n(4) &
                  & + k(164)*n(8)

      Jmtx_tmp(7,15) = k(111)*n(8) &
                  & + k(119)*n(17)

      Jmtx_tmp(7,16) = k(127)*n(8)

      Jmtx_tmp(7,19) = k(157)*n(8)

      Jmtx_tmp(8,1) = k(199)*n(19) &
                  & - k(192)*n(14)

      Jmtx_tmp(8,2) = k(97)*n(17)

      Jmtx_tmp(8,3) = k(71)*n(17)

      Jmtx_tmp(8,4) = k(86)*n(13) &
                  & + k(91)*n(17)

      Jmtx_tmp(8,5) = k(146)*n(13) &
                  & + k(149)*n(17)

      Jmtx_tmp(8,7) = k(175)*n(13) &
                  & + k(178)*n(17) &
                  & + k(217)*n(4)*n(6)

      Jmtx_tmp(8,8) = - k(39)*n(3) &
                  & - k(40)*n(3) &
                  & - k(192)*n(1) &
                  & - k(193)*n(11) &
                  & - k(194)*n(11) &
                  & - k(195)*n(13) &
                  & - k(196)*n(13) &
                  & - k(197)*n(17) &
                  & - k(198)*n(17)

      Jmtx_tmp(8,11) = k(209)*n(13)

      Jmtx_tmp(8,12) = k(199)*n(1)

      Jmtx_tmp(8,13) = k(102)*n(13)

      Jmtx_tmp(8,15) = k(115)*n(13)

      Jmtx_tmp(8,16) = k(125)*n(8)

      Jmtx_tmp(8,17) = k(131)*n(8)

      Jmtx_tmp(9,1) = k(192)*n(14) &
                  & - k(185)*n(15)

      Jmtx_tmp(9,2) = k(96)*n(17)

      Jmtx_tmp(9,3) = k(67)*n(11) &
                  & + k(72)*n(17)

      Jmtx_tmp(9,4) = k(87)*n(13) &
                  & + k(92)*n(17)

      Jmtx_tmp(9,5) = k(144)*n(11)

      Jmtx_tmp(9,6) = k(135)*n(13)

      Jmtx_tmp(9,7) = k(170)*n(4)

      Jmtx_tmp(9,8) = k(192)*n(1)

      Jmtx_tmp(9,9) = - k(37)*n(3) &
                  & - k(38)*n(3) &
                  & - k(185)*n(1) &
                  & - k(186)*n(8) &
                  & - k(187)*n(11) &
                  & - k(188)*n(13) &
                  & - k(189)*n(17) &
                  & - k(190)*n(17) &
                  & - k(191)*n(17)

      Jmtx_tmp(9,13) = k(103)*n(13) &
                  & + k(104)*n(17)

      Jmtx_tmp(9,14) = k(107)*n(11) &
                  & + k(109)*n(13)

      Jmtx_tmp(9,15) = k(112)*n(8) &
                  & + k(114)*n(13) &
                  & + k(120)*n(17)

      Jmtx_tmp(9,16) = k(126)*n(8)

      Jmtx_tmp(9,18) = k(151)*n(11)

      Jmtx_tmp(9,20) = k(212)*n(11)

      Jmtx_tmp(10,4) = k(83)*n(11) &
                  & + k(89)*n(13)

      Jmtx_tmp(10,10) = - k(33)*n(3) &
                  & - k(34)*n(3) &
                  & - k(163)*n(4) &
                  & - k(164)*n(8) &
                  & - k(165)*n(8) &
                  & - k(166)*n(8) &
                  & - k(167)*n(8) &
                  & - k(168)*n(11) &
                  & - k(169)*n(13)

      Jmtx_tmp(10,19) = k(155)*n(4) &
                  & + k(156)*n(8)

      Jmtx_tmp(11,1) = - k(205)*n(18)

      Jmtx_tmp(11,3) = k(69)*n(17)

      Jmtx_tmp(11,11) = - k(43)*n(3) &
                  & - k(44)*n(3) &
                  & - k(205)*n(1) &
                  & - k(206)*n(11) &
                  & - k(207)*n(11) &
                  & - k(208)*n(11) &
                  & - k(209)*n(13) &
                  & - k(210)*n(17) &
                  & - k(211)*n(17)

      Jmtx_tmp(12,1) = k(205)*n(18) &
                  & - k(199)*n(19)

      Jmtx_tmp(12,2) = k(98)*n(17)

      Jmtx_tmp(12,3) = k(70)*n(17)

      Jmtx_tmp(12,5) = k(147)*n(13)

      Jmtx_tmp(12,6) = k(133)*n(8) &
                  & + k(138)*n(17)

      Jmtx_tmp(12,7) = k(179)*n(17)

      Jmtx_tmp(12,9) = k(188)*n(13) &
                  & + k(189)*n(17)

      Jmtx_tmp(12,11) = k(205)*n(1) &
                  & + k(206)*n(11)

      Jmtx_tmp(12,12) = - k(41)*n(3) &
                  & - k(42)*n(3) &
                  & - k(199)*n(1) &
                  & - k(200)*n(8) &
                  & - k(201)*n(11) &
                  & - k(202)*n(11) &
                  & - k(203)*n(13) &
                  & - k(204)*n(17)

      Jmtx_tmp(12,13) = k(105)*n(17)

      Jmtx_tmp(12,14) = k(108)*n(13) &
                  & + k(110)*n(17)

      Jmtx_tmp(12,15) = k(121)*n(17)

      Jmtx_tmp(12,17) = k(130)*n(8)

      Jmtx_tmp(12,18) = k(152)*n(13) &
                  & + k(153)*n(17)

      Jmtx_tmp(12,20) = k(213)*n(13)

      Jmtx_tmp(13,1) = - k(101)*n(20)

      Jmtx_tmp(13,3) = k(63)*n(6)

      Jmtx_tmp(13,13) = - k(20)*n(3) &
                  & - k(100)*n(4) &
                  & - k(101)*n(1) &
                  & - k(102)*n(13) &
                  & - k(103)*n(13) &
                  & - k(104)*n(17) &
                  & - k(105)*n(17)

      Jmtx_tmp(14,2) = k(214)*n(4)*n(4)

      Jmtx_tmp(14,3) = k(61)*n(4)

      Jmtx_tmp(14,13) = k(100)*n(4)

      Jmtx_tmp(14,14) = - k(22)*n(3) &
                  & - k(23)*n(3) &
                  & - k(106)*n(8) &
                  & - k(107)*n(11) &
                  & - k(108)*n(13) &
                  & - k(109)*n(13) &
                  & - k(110)*n(17)

      Jmtx_tmp(15,1) = k(123)*n(24)

      Jmtx_tmp(15,15) = - k(24)*n(3) &
                  & - k(111)*n(8) &
                  & - k(112)*n(8) &
                  & - k(113)*n(11) &
                  & - k(114)*n(13) &
                  & - k(115)*n(13) &
                  & - k(116)*n(13) &
                  & - k(117)*n(13) &
                  & - k(118)*n(13) &
                  & - k(119)*n(17) &
                  & - k(120)*n(17) &
                  & - k(121)*n(17) &
                  & - k(122)*n(17)

      Jmtx_tmp(15,16) = k(123)*n(1)

      Jmtx_tmp(16,1) = - k(123)*n(24)

      Jmtx_tmp(16,4) = k(78)*n(8) &
                  & + k(85)*n(11)

      Jmtx_tmp(16,16) = - k(25)*n(3) &
                  & - k(123)*n(1) &
                  & - k(124)*n(4) &
                  & - k(125)*n(8) &
                  & - k(126)*n(8) &
                  & - k(127)*n(8) &
                  & - k(128)*n(11)

      Jmtx_tmp(17,4) = k(79)*n(8) &
                  & + k(90)*n(13)

      Jmtx_tmp(17,16) = k(124)*n(4)

      Jmtx_tmp(17,17) = - k(26)*n(3) &
                  & - k(129)*n(4) &
                  & - k(130)*n(8) &
                  & - k(131)*n(8) &
                  & - k(132)*n(11)

      Jmtx_tmp(18,1) = - k(150)*n(29)

      Jmtx_tmp(18,3) = k(64)*n(8)

      Jmtx_tmp(18,5) = k(141)*n(4) &
                  & + k(142)*n(8)

      Jmtx_tmp(18,6) = k(215)*n(4)*n(4) &
                  & + k(216)*n(4)*n(6)

      Jmtx_tmp(18,14) = k(106)*n(8)

      Jmtx_tmp(18,18) = - k(30)*n(3) &
                  & - k(31)*n(3) &
                  & - k(150)*n(1) &
                  & - k(151)*n(11) &
                  & - k(152)*n(13) &
                  & - k(153)*n(17) &
                  & - k(154)*n(17)

      Jmtx_tmp(19,4) = k(84)*n(11)

      Jmtx_tmp(19,19) = - k(32)*n(3) &
                  & - k(155)*n(4) &
                  & - k(156)*n(8) &
                  & - k(157)*n(8) &
                  & - k(158)*n(8) &
                  & - k(159)*n(8) &
                  & - k(160)*n(8) &
                  & - k(161)*n(11) &
                  & - k(162)*n(13)

      Jmtx_tmp(20,18) = k(154)*n(17)

      Jmtx_tmp(20,20) = - k(45)*n(3) &
                  & - k(212)*n(11) &
                  & - k(213)*n(13)

      Jmtx_tmp(21,15) = k(113)*n(11) &
                  & + k(116)*n(13)

      Jmtx_tmp(21,19) = k(158)*n(8)

      Jmtx_tmp(21,21) = - k(46)*n(3)

      Jmtx_tmp(22,15) = k(117)*n(13)

      Jmtx_tmp(22,16) = k(128)*n(11)

      Jmtx_tmp(22,19) = k(159)*n(8)

      Jmtx_tmp(22,22) = - k(47)*n(3)

      Jmtx_tmp(23,5) = k(145)*n(11)

      Jmtx_tmp(23,6) = k(134)*n(11) &
                  & + k(136)*n(13)

      Jmtx_tmp(23,7) = k(176)*n(13) &
                  & + k(180)*n(17)

      Jmtx_tmp(23,8) = k(193)*n(11)

      Jmtx_tmp(23,10) = k(165)*n(8)

      Jmtx_tmp(23,12) = k(201)*n(11)

      Jmtx_tmp(23,15) = k(118)*n(13) &
                  & + k(122)*n(17)

      Jmtx_tmp(23,17) = k(132)*n(11)

      Jmtx_tmp(23,19) = k(160)*n(8)

      Jmtx_tmp(23,23) = - k(48)*n(3)

      Jmtx_tmp(24,7) = k(171)*n(8) &
                  & + k(181)*n(17)

      Jmtx_tmp(24,10) = k(166)*n(8)

      Jmtx_tmp(24,24) = - k(49)*n(3)

      Jmtx_tmp(25,5) = k(148)*n(13)

      Jmtx_tmp(25,6) = k(137)*n(13) &
                  & + k(139)*n(17)

      Jmtx_tmp(25,7) = k(172)*n(8) &
                  & + k(182)*n(17)

      Jmtx_tmp(25,8) = k(195)*n(13)

      Jmtx_tmp(25,9) = k(186)*n(8) &
                  & + k(190)*n(17)

      Jmtx_tmp(25,10) = k(167)*n(8)

      Jmtx_tmp(25,11) = k(207)*n(11)

      Jmtx_tmp(25,12) = k(203)*n(13)

      Jmtx_tmp(25,25) = - k(50)*n(3)

      Jmtx_tmp(26,8) = k(197)*n(17)

      Jmtx_tmp(26,26) = - k(51)*n(3)

      Jmtx_tmp(27,6) = k(140)*n(17)

      Jmtx_tmp(27,8) = k(198)*n(17)

      Jmtx_tmp(27,12) = k(200)*n(8)

      Jmtx_tmp(27,27) = - k(52)*n(3)

      Jmtx_tmp(28,11) = k(210)*n(17)

      Jmtx_tmp(28,28) = - k(53)*n(3)

      Jmtx_tmp(29,11) = k(211)*n(17)

      Jmtx_tmp(29,29) = - k(54)*n(3)

      Jmtx_tmp(30,19) = k(161)*n(11)

      Jmtx_tmp(30,30) = - k(55)*n(3)

      Jmtx_tmp(31,7) = k(173)*n(11)

      Jmtx_tmp(31,10) = k(168)*n(11)

      Jmtx_tmp(31,31) = - k(56)*n(3)

      Jmtx_tmp(32,7) = k(174)*n(11)

      Jmtx_tmp(32,9) = k(187)*n(11)

      Jmtx_tmp(32,32) = - k(57)*n(3)

      Jmtx_tmp(33,7) = k(177)*n(13) &
                  & + k(183)*n(17)

      Jmtx_tmp(33,8) = k(194)*n(11)

      Jmtx_tmp(33,12) = k(202)*n(11)

      Jmtx_tmp(33,33) = - k(58)*n(3)

      Jmtx_tmp(34,7) = k(184)*n(17)

      Jmtx_tmp(34,8) = k(196)*n(13)

      Jmtx_tmp(34,9) = k(191)*n(17)

      Jmtx_tmp(34,11) = k(208)*n(11)

      Jmtx_tmp(34,34) = - k(59)*n(3)

      Jmtx_tmp(35,12) = k(204)*n(17)

      Jmtx_tmp(35,35) = - k(60)*n(3)


      ! Converting index --------------------------
      if (nz >= 2) then
        do i1 = 1, nsp_var_out
          jsp = var_out_all_in(i1)
          i0 = all_in_var_in(jsp)
          do j1 = 1, nsp_var_out
            jsp = var_out_all_in(j1)
            j0 = all_in_var_in(jsp)
            i = (iz-1)*nsp_var_out+i1
            j = (iz-1)*nsp_var_out+j1 + nsp_var_out + 1 - i
            Jmtx(j,i) = - Jmtx_tmp(i0,j0)
          end do 
        end do 
      else if (nz == 1) then
        do i1 = 1, nsp_var_out
          jsp = var_out_all_in(i1)
          i0 = all_in_var_in(jsp)
          do j1 = 1, nsp_var_out
            jsp = var_out_all_in(j1)
            j0 = all_in_var_in(jsp)
            Jmtx(j0,i0) = - Jmtx_tmp(i0,j0)
          end do 
        end do 
      end if

    end do 

  end subroutine p__PROTEUS_Jacobian


  subroutine get_solar_flux_data(lambda, dlambda, nwl, fname, unit1, unit2, & ! in
    &                            solar_flux                                 ) ! out
    implicit none
    real(dp), intent(in)  :: lambda(1:), dlambda(1:)
    integer,  intent(in)  :: nwl
    character(len=*), intent(in) :: unit1, unit2, fname
    real(dp), intent(out) :: solar_flux(1:)
    integer i, nh, il, nl
    real(dp), allocatable :: idata(:,:), odata(:)

    write(*,'(a)',advance='no')  '  Reading datafile: '//trim(ADJUSTL(fname))//'...'

    allocate(odata(nwl))

    call get_header_line_number(nh, fname)
    call get_data_line_number(nl, fname)
    allocate(idata(nl,2))

    open(11, file = fname, status = 'old')
      do i = 1, nh; read(11,*); end do
      do il = 1, nl
        read(11,*) idata(il,1), idata(il,2) 

        ! Unit conversion: column 1
        if (unit1/='nm' .and. unit1/='A' .and. unit1/='/cm' .and. unit1/='cm^-1') then 
          print *, ''
          print *, 'error!'
          print *, '  The unit of the datafile "'//trim(adjustl(fname))//'" column 1 is not recognized.'
          print *, '  The unit of column 1 should be "nm", "A", "/cm" or "cm^-1".'
          stop
        end if
        if (unit1 == 'A') then 
          idata(il,1) = idata(il,1) * 1.0e-1_dp ! [A -> nm]
        else if (unit1 == '/cm' .or. unit1 == 'cm^-1') then 
          idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
        end if

        ! Unit conversion: column 2
        if (unit2 == '/cm^2/s/nm') then 
          idata(il,2) = idata(il,2) * 1.0e4_dp ! [/cm^2/s/nm -> /m^2/s/nm]
        else if (unit2 == 'W/m^2/nm') then 
          idata(il,2) = idata(il,2) * idata(il,1) * 1.0e-9_dp / 6.626e-34_dp / 2.99792458e8_dp ! [W/m^2/nm  ->  /m^2/s/nm]
        end if

      end do
    close(11)

    write(*,*) 'done.'

    call binning_flux(lambda, dlambda, idata, nl, nwl, odata)
    solar_flux(1:nwl) = odata(1:nwl)
    deallocate(idata)
    deallocate(odata)

  end subroutine get_solar_flux_data


  subroutine binning_flux(lambda, dlambda, idata, nl, nwl, & ! in
    &                     odata                            ) ! out
    implicit none
    real(dp),   intent(in)  :: lambda(nwl), dlambda(nwl)
    integer,    intent(in)  :: nl, nwl
    real(dp),   intent(in)  :: idata(nl,2)
    real(dp),   intent(out) :: odata(nwl)
    integer iwl, il, label, il0, ndata
    real(dp) idata_tmp(nl,2)
    real(dp) i0, ip, im, dip, dim, o0, dop, dom, sum_data, sum_wl, tmp
  
    if (idata(1,1) > idata(nl,1)) then 
      do il = 1, nl
        idata_tmp(il,1) = idata(nl-il+1,1)
        idata_tmp(il,2) = idata(nl-il+1,2)
      end do 
    else 
      idata_tmp(:,:) = idata(:,:)
    end if
  
    odata = 0.0_dp
    il0 = 1
  
    do iwl = 1, nwl 
  
      label = -1
      ndata = 0
      sum_data = 0.0_dp
      sum_wl   = 0.0_dp
      tmp = 0.0_dp
  
      loop: do il = il0, nl 
        
        ! input bin
        if (il == 1) then 
          i0 = idata_tmp(il,1)
          ip = idata_tmp(il+1,1)
          dip = ip - i0
          dim = dip
        else if (il > 1 .and. il < nl) then 
          im = idata_tmp(il-1,1)
          i0 = idata_tmp(il,1)
          ip = idata_tmp(il+1,1)
          dim = i0 - im
          dip = ip - i0
        else if (il == nl) then 
          im = idata_tmp(il-1,1)
          i0 = idata_tmp(il,1)
          dim = i0 - im
          dip = dim
        end if
  
        ! model bin
        if (iwl == 1) then 
          o0 = lambda(1)
          dop = dlambda(1) * 0.5d0
          dom = dlambda(1) * 0.5d0
        else if (iwl > 1 .and. iwl < nwl) then 
          o0 = lambda(iwl)
          dom = dlambda(iwl) * 0.5d0
          dop = dlambda(iwl) * 0.5d0
        else if (iwl == nwl) then 
          o0 = lambda(iwl)
          dom = dlambda(iwl) * 0.5d0
          dop = 1.0e10_dp ! to take large enough value
        end if 
        
        ! for binning
        if (o0-dom <= i0 .and. i0 <= o0+dop) then 
          sum_data = sum_data + idata_tmp(il,2)*(dim+dip)
          sum_wl   = sum_wl + (dim+dip)
          ndata = ndata+1
          if (ndata >= 2) label = 1
          if (iwl==nwl .or. il==nl) label = 1
          if (iwl==1 .or. il==1) label = 1
        end if
        
        ! for interpolate
        if (i0 <= o0 .and. o0 < ip) then 
          tmp = (idata_tmp(il,2)*(ip-o0) + idata_tmp(il+1,2)*(o0-i0))/(ip-i0)
          if (label /= 1) label = 0
        end if
  
        ! to escape from this loop
        if (label == 1 .and. i0 > o0+dop) then 
          if (il >= 2) il0 = il-1
          exit loop
        else if (label == 0 .and. ip < o0) then 
          if (il >= 2) il0 = il-1
          exit loop
        end if
        
      end do loop ! il
  
      ! binning if an input bin is smaller than a model bin
      if (label == 1) then 
        odata(iwl) = sum_data / sum_wl
      end if
  
      ! interpolate if an input bin is larger than a model bin
      if (label == 0) then 
        odata(iwl) = tmp
      end if
  
    end do ! iwl
  
  end subroutine binning_flux
  
  
  subroutine get_header_line_number(nh, fname)
    implicit none
    integer nh
    character(len=256) fname
    character(len=256) strm
    open(11, file = fname, status = 'old' )
      nh = 0
      do
        read(11,'(A)') strm
        if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
          nh = nh + 1
        else
          exit
        end if
      end do
    close(11)
  end subroutine get_header_line_number
  
  
  subroutine get_data_line_number(nl, fname)
    implicit none
    integer nl, ios
    character(len=256) fname
    character(len=256) strm
    open(11, file = fname, status = 'old' )
      nl = 0
      do
        read(11,'(A)',iostat = ios) strm
        if (ios < 0) exit
        if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
          cycle
        else
          nl = nl + 1
        endif
      end do
    close(11)
  end subroutine get_data_line_number
  
  
end module p__PROTEUS
