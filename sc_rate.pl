#!/usr/bin/perl
# Last Update by ./update_subroutines.pl: Thu Jan  2 03:26:56 GMT 1997
#________________________________________________________________________
# Title     : sc_rate.pl  (or sc_rate in certain cases)
# Usage     : sc_rate.pl test_file.msf [test_structure_file.jp]
# Function  : Evaluates the quality of any sequence alignment.
#             This caculates, basically, the Sequence identity(in %) divided
#             by Composition identity rate of sequences with certain scanning
#             window size of (8 ~ 24)
# Example   : sc_rate.pl aa.msf aa.jp
# Warning   : Takes 'ALIGNED' sequences. So, need at least a pair of seqs aligned.
# Class     : Bio
# Keywords  : SC rate, Sequence identity, composition identity
# Options   : Following vars are read into the program at initialization by &parse_arguments
#             subroutine, so please leave them. It is not just comments.
#    $av_sc_segment     becomes    1  by  a  # To smooth the SC rates
#                                  SC rates with Window size defined by a=X  .
#    $make_gap_dot      becomes    d  by  d  # Modifies SC rate by a dynamic factor of window
#    $av_sc_segment           =   [ ] by  a= # window size for averaging after getting SC rates
#    $show_calculation  becomes    s  by  s  # Shows the scanning windows during calculation
#    $no_stat           becomes    ns by  ns # do not show stat with 2 file inputs(
#    $variable_win_size becomes    v  by  v  # This sets the scan win size dynamic
#    $output_ordered    becomes    o  by  o  # in final printing, seq
#    $interlaced_out    becomes    i  by  i  # in final printing, seq
#    $normal_sc_rate    becomes    n  by  n  # Normalize (0-9) the SC rates
#    $normal_er_rate    becomes    e  by  e  # Normalize (0-9) the error rates(2 file input)
#    $conv_err_2_bin    becomes    c  by  c  # Converts error rates to 0 or 1 (occurred or not)
#    $redu_window       becomes    r  by  r  # Tries to reduce the scan wind
#                                            This is only activated with 'v' option.
#    $apply_factor      becomes    f  by  f  # Modifies SC rate by a dynamic factor of window
#                                            size. This is only activated with 'v' option
#    $minus_whole_cs    becomes    m  by  m  # substracts rates with whole seq
#                                            # this makes rates more strict in predicting
#    $print_width_size  is        [ ] by  w= # the final output print block width
#
#    $ss_opt            becomes    ss by  ss, SS, -ss, -SS     #  for secondary structure only
#    $H                 =         'H' by   -H or -h or H or h  # to retrieve only H segment
#    $S                 becomes   'S' by   -S or  S            # to retrieve only S segment
#    $E                 becomes   'E' by   -E or  E            # to retrieve only E segment
#    $T                 becomes   'T' by   -T or -t or T or t  # to retrieve only T segment
#    $I                 becomes   'I' by   -I or  I            # to retrieve only I segment
#    $G                 becomes   'G' by   -G or -g or G or g  # to retrieve only G segment
#    $B                 becomes   'B' by   -B or -b or B or b  # to retrieve only B segment
#    $HELP              becomes    1  by   -help   # for showing help
#    $simplify          becomes    1  by   -p or P or -P, p
#    $simplify          becomes    1  by   -simplify or simplify, Simplify SIMPLIFY
#    $comm_col          becomes   'C' by   -C or C or common
#    $LIMIT             becomes    L  by   -L, L               # to limit the error rate to 9 .
#    $HELP              becomes    1  by   -h, h, '?'          # for showing help
# Package   :
# Reference :
# Returns   : No file writing.
# Tips      : This is handy to see how reliable your alignment is.
# Argument  : Upto 2 files and various options (in '-h'  or 'h' format )
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   : Bio
# Enclosed  :
#--------------------------------------------------------------------

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#     INITIALIZATION and PARSING ARGUMENTS AND OPTIONS
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

my(@input_files)=@{&parse_arguments(1)};
if( $debug eq 1){  # to activate debugging, you just put '#' or '_'(for bash) at prompt
  print "\n",__LINE__, " Your input file(s) are :  ", "@input_files\n\n" if $debug eq 1;
  print "\n",__LINE__, " \$ss_opt is set to :  ", "$ss_opt\n";
  print "\n",__LINE__, " \$av_sc_segment is set to :  ", "$av_sc_segment\n";
  print "\n",__LINE__, " \$make_gap_dot is set to :  ", "$make_gap_dot\n";
}
	#""""""""""""""""""""""""""
	#        M A I N
	#__________________________

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# When there is only one input seq. alignment file (i.e. just computer alignment is given)
#_______________________________________________________________________________________
if( @input_files == 1 ){
   $base =${&get_base_names(\$input_files[0])};
   my(%input)= %{&read_any_seq_files(\$input_files[0])};


   ############# scan_win_and_get_sc_rate_pairs is MAIN routine ###############
   my(%out1) = %{&scan_win_and_get_sc_rate_pairs(\%input, \$window_size,
			   \$variable_win_size, \$show_calculation, \$minus_whole_cs,
			   \$redu_window, \$apply_factor, \$make_gap_dot)};
   ############# scan_win_and_get_sc_rate_pairs is MAIN routine ###############


   &print_seq_in_block(\%out1,  \%input, \$print_width_size);
	  print "\n To know how to interpret it, press 'i' and 'return' or just 'return'\n";
   $key_in = getc;
   if($key_in eq 'i'){
	  print "\n\n With  xxxxxx(2), the (2) means the ratio of Compos. id divided by Seq. id. \n";
	  print " Therefore, the higher the worse the alignment will be.\n";
	  print "\n With 33333333332122333444354344453554332222 , you can tell how reliable the";
	  print "\n alignment at the residue position, 9 => best reliability, 0 => worst\n\n\n";
   }else{
	 print "\n Bye \n\n\n";
   }
}elsif( @input_files > 1 ){
	print "\n",__LINE__, " 2 input file, getting error rate too.  \n" if $debug eq 1;
	$base =${&get_base_names(\$input_files[0])};

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	##               VAR   values initialization
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	$av_sc_segment = 4 if !defined($av_sc_segment);  # smoothing window size
	$threshold = 1; # Threshold '1' means everything above 0 is 1 in error rates.

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	##             Getting input file order right                             ##
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if($input_files[0]=~/\.msf/){  $input = $input_files[0];  $input2 = $input_files[1]; }
	elsif($input_files[0]=~/\.jp/){ $input = $input_files[1]; $input2 = $input_files[0]; }

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#             Reading input files                                        ##
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my(%input)= %{&read_any_seq_files(\$input)};
	my(%input2)=%{&read_any_seq_files(\$input2)};

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#             Getting error rates and SC rates                          ###
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my(%error_rate)=%{&get_position_shift_rate(\%input, \%input2, $ss_opt)};

	print "\n",__LINE__," $0 : \$ss_opt was passed as \"$ss_opt\"\n" if $debug eq 1;
	print "\n",__LINE__," $0 : This produces error rates for sec. str. regions only\n" if $debug eq 1;
	print "\n",__LINE__," The error rate hash is\n", &show_hash(\%error_rate) if $debug eq 1;
	my(%sc_rate)  = %{&scan_win_and_get_sc_rate_pairs(\%input, \$window_size,
	  \$make_gap_dot, \$variable_win_size, \$show_calculation, \$redu_window,
	  \$apply_factor, \$minus_whole_cs)};

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#              Sets the option inputs vars                                ##
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if( $conv_err_2_bin =~/[cC]/ ){
	   %error_rate=%{ &convert_num_0_or_1_hash_opposite(\%error_rate, \$threshold)};
	}
	if( $normal_er_rate =~ /[eE]+/ ){  ## option given by 'e'
	   %error_rate=%{ &normalize_numbers(\%error_rate)};
	}
	if($normal_sc_rate =~ /[nN]+/ ){   ## option given by 'n'
	   %sc_rate =  %{&normalize_numbers(\%sc_rate)};
	}
	if( $av_sc_segment =~ /^\d+/ ){
	   %sc_rate  = %{&scan_win_get_average(\%sc_rate, $av_sc_segment)};
	}

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#     Getting the statistics of the reliability compared to error rate    ##
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	unless( $no_stat =~ /^ns/ ){   #### if the 'ns' option is not set.
	   ($tally_result, $num_occurances)=&tally_2_hashes(\%sc_rate, \%error_rate); # n is for NOT normalize
	   %tally_result=%{$tally_result};
	   %tally_num_occurances=%{$num_occurances}; ## This is to show how many '0','1', etc has occurred.
	   print "\n";

	   ($tally_result, $num_occurances)=&tally_2_hashes(\%sc_rate, \%error_rate, 'n', 'i'); # n is for NOT normalize
	   %tally_result_average=%{$tally_result}; # n is for NOT normalize
	   ($tally_result, $num_occurances)=&tally_2_hashes(\%sc_rate, \%error_rate, 'a'); # n is for NOT normalize
	   %tally_result_percent=%{$tally_result}; # n is for NOT normalize
	}

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#    Final  Printing OUT to SCREEN
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	print "\n";
	&print_seq_in_block(\%input, \$print_width_size);
	#%sc_rate =%{&superpose_seq_hash(\%sc_rate, \%error_rate)};
	&print_seq_in_block(\%sc_rate, \%error_rate, \$interlaced_out, \$print_width_size);  #, \%error_rate  \%sst_hash,
	%header1=  ('SC', 'AddedUp ');
	%header1_1=('SC', 'Occur   ');
	%header2=  ('SC', 'Err Av  ');
	%header3=  ('SC', 'Ratio   ');
	&print_seq_in_columns(\%header1, \%header1_1, \%header2,\%header3 );
	&print_seq_in_columns(\%tally_result, \%tally_num_occurances,
						  \%tally_result_percent, \%tally_result_average, 's11');
}

		#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		####  Sub Routines
		#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



#________________________________________________________________________
# Title     : scan_win_and_get_sc_rate_pairs
# Usage     : %out1 = %{&scan_win_and_get_sc_rate_pairs(\%input, \$window_size)};
# Function  : scans input sequences(arg1) in a given(arg2) window size and gets
#             each composition and sequence identity rate(sc_rate) of the window.
#             sc rate = Sequence Id(%)/ Composition Id(%)
# Example   : input hash: ( seq1,  'ABCDEFG.HIK',      (2 or more sequences accepted)
#                         seq2,  'DFD..ASDFAFS',
#                         seq3,  'DDDDD..ASDFAFS' );
#             input winsize : 5;
#
#             output hash; (seq1seq2, 1,2,2,2,1,1,2,2); <-- joined by ',';
#             output hash; (seq1seq3, 1,2,2,2,1,1,2,2); <-- joined by ',';
#                  The numbers are ratios(compos/seqid) with given window size.
# Warning   : when $seqid is zero  the rate becomes $compos_id/10   !!!
# Keywords  :
# Options   :
# Returns   : a reference of a hash.
# Argument  : One ref. for hash, one ref. for a scalar.
# Version   : 1.0
#--------------------------------------------------------------------
sub scan_win_and_get_sc_rate_pairs{
  my($base_l)=1; my($scale)=1;
  my(%input, $i, $j, $window_size, $show_calculation, $redu_window,
	  $variable_win_size);
  for($i=0; $i < @_; $i ++){
	 if( (ref($_[$i])) eq 'HASH'  ){ %input=%{$_[$i]};}
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /^d+/)){ $window_size= ${$_[$i]}; }
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /^[vV]+/)){ $variable_win_size= 'v' }
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /[sS]+/)){ $show_calculation = 's'}
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /[rR]+/)){ $redu_window  = 'r'}
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /[fF]+/)){ $apply_factor  = 'f'}
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /[dD]+/)){ $make_gap_dot  = 'd'}
	 elsif( (ref($_[$i]) eq 'SCALAR' )&&(${$_[$i]} =~ /[mM]+/)){ $minus_whole_cs  = 'm'}
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^\d+/) ){ $window_size= $_[$i]; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[vV]+/) ){ $variable_win_size = 'v'; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[sS]+/) ){ $show_calculation  = 's'; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[rR]+/) ){ $redu_window  = 'r'; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[fF]+/) ){ $apply_factor  = 'f'; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[dD]+/) ){ $make_gap_dot  = 'd'; }
	 elsif( (!(ref($_[$i])))&&($_[$i] =~ /^[mM]+/) ){ $minus_whole_cs  = 'm'; }
  }

  if(defined(${$_[2]})){ $base_l=${$_[2]}; }
  if(defined(${$_[3]})){ $scale =${$_[3]}; }
  my(@sequences, @out_rate, $i, $j, $title, $window_1, $window_2,
	  $ratio_compos_vs_seqid, @array_of_2_seq,%out_hash );
  my(@keys)= sort keys %input;
  for ($i=0; $i<@keys; $i++){    # putting sequences from hash to an array
	 for($j=$i+1; $j< @keys; $j++){
		push(@sequences, $input{$keys[$i]}, $input{$keys[$j]});

		######################################################################
		#################  PASSING OVER TO the next SUB routine ##############
		######################################################################
		#---> @sequences will have ('ABCDEFG.HIK', 'DFD..ASDFAFS'); ##########
		($out_rate_arr_ref,$whole_rate_ref)= &get_windows_sc_rate_array(

						\@sequences,\$window_size, $variable_win_size, \$apply_factor,
						\$redu_window, $make_gap_dot, $show_calculation, \$minus_whole_cs );

		undef(@sequences);
		@out_rate=@{$out_rate_arr_ref};
		$whole_rate=${$whole_rate_ref};
		$title = "$keys[$i]_$keys[$j]\($whole_rate\)";
		$out_hash{$title}=join(",", @out_rate);
	 }
  }
  return( \%out_hash );
}

#________________________________________________________________________
# Title     : get_windows_sc_rate_array
# Usage     : @out_rate = @{&get_windows_compos_and_seqid_rate_array(\@seq, \$win_size)};
# Function  : actual working part of scan_windows_and_get_compos_seqid_rate
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : \@ratio_array, \$ratio_whole_seq
# Tips      :
# Argument  : (\@input, \$window_size);  @input => ('ABCDEFG.HIK', 'DFD..ASDFAFS', 'ASDFASDFASAS');
#             Input ar => ( 'ABCDEFG
#                'DFD..ASDFAFS'
#                'ASDFASDFASAS' )  as the name of  @sequences.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_windows_sc_rate_array{
  my($base_level)=1; my($scale)=1; my($window_size, $show_calculation, $redu_window,
	 @input, @input0, @input1, $variable_win_size, $apply_factor,
	 @ratio_array, @array_of_2_seq, $seq_id, $offset, $half_of_w_size, $t, $length, $w,
	 $compos_id, $seq_id, $window_2, $window_1, $compos_whole_seq, $seq_id_whole_seq,
	 $ratio_whole_seq, $win_rate_div_by_whole_rate, $normalizing_factor, $lowest_rate,
	 $winsize_reduc_factor, $largest_win_reached, $ori_win_size);

	#""""""""""""""""""""""< handle_arguments{ head Ver 1.3 >"""""""""""""""""""
	my(@A )=&handle_arguments( @_ );my( $num_opt )=${$A[7]};my( $char_opt )=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};my(@raw_string)=@{$A[9]};
	my($i, $j, $c, $d, $e, $f, $g, $h, $k, $l, $p, $q, $r, $s, $t, $u, $v, $w, $x,$y,$z);
	if($debug==1){ print "\n   \@hash has \"@hash\"\n   \@raw_string has   \"@raw_string\"
	\@array has \"@array\"\n   \@char_opt has   \"@char_opt\"\n   \@file has \"@file\"
	\@string has \"@string\""; }
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   if( $num_opt =~/^(d+)/){  $window_size= $1; }
   if( $char_opt=~/v/i){  $variable_win_size= 'v' }
   if( $char_opt=~/s/i){  $show_calculation = 's'}
   if( $char_opt=~/r/i){  $redu_window  = 'r'}
   if( $char_opt=~/f/i){  $apply_factor  = 'f'}
   if( $char_opt=~/d/i){  $make_gap_dot  = 'd'}
   if( $char_opt=~/m/i){  $minus_whole_cs  = 'm'}

  @input = @{$array[0]};
  if(defined(${$_[2]})){ $base_level =${$_[2]}; }
  if(defined(${$_[3]})){ $scale  =${$_[3]}; }

  for ($t=0; $t< @input; $t++){ $length=length($input[$t]) if(length($input[$t])>$length);}
  if ($length < $window_size){  $window_size = $length;   }
	 #___________ getting ratio for the whole sequence ___________

  $compos_whole_seq=${&main::compos_id_percent_array(\@input)};  ## for whole composition rate
  $seq_id_whole_seq=${&main::seq_id_percent_array(\@input)};
  print "\nComposition ID of the alignment:  $compos_whole_seq\%\n";
  print "Sequence    ID of the alignment:  $seq_id_whole_seq\%\n";
  if ($seq_id_whole_seq == 0){  $ratio_whole_seq =0; }
  else{  $ratio_whole_seq = $compos_whole_seq/$seq_id_whole_seq;  }
  print "Composition and Sequ.  ID Ratio:   $ratio_whole_seq\n";

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ##########   Initial Window size setting      ##############################
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if   (  $ratio_whole_seq >9 ){ $window_size= 28; } # $ratio_whole_seq is
  elsif(  $ratio_whole_seq >8 ){ $window_size= 26; } # a whole CS rate !!
  elsif(  $ratio_whole_seq >7 ){ $window_size= 24; }
  elsif(  $ratio_whole_seq >6 ){ $window_size= 22; }
  elsif(  $ratio_whole_seq >5 ){ $window_size= 20; }
  elsif(  $ratio_whole_seq >4 ){ $window_size= 18; }
  elsif(  $ratio_whole_seq >3 ){ $window_size= 16; }
  elsif(  $ratio_whole_seq >2 ){ $window_size= 12; }
  elsif(  $ratio_whole_seq >0 ){ $window_size= 8; }

  print "Window size used is :  $window_size\n";
  #$window_size = 10;
  $largest_win_reached = 24;
  $ori_win_size        = $window_size;

  #----------- Spliting the seq. into arrays to enable $make_gap_dot var -----
  if( $make_gap_dot =~ /^[dD]+/){
	 $input[0] =~s/,//g;             $input[1] =~s/,//g;
	 @input0 = split('', $input[0]); @input1 = split('', $input[1]);
  }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #        MAIN Calc part
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for ($w=0; $w < $length; $w++){

	 $largest_win_reached = $window_size if $window_size > $largest_win_reached;

	 #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###          FACTOR calculation                                    ##
	 #####################################################################
	 $factor= ($window_size/40)*2  if ($apply_factor eq 'f');

	 #####################################################################
	 ##       Turn RATES to DOTs when there are gaps                    ##
	 #####################################################################
	 # @input0, @input1 are the whole length sequence splited.
	 if($make_gap_dot =~ /^[dD]+/){
		if( ($input0[$w] eq '.') || ($input1[$w] eq '.') ){
		  $ratio_compos_vs_seqid = '.';
		  push(@ratio_array, $ratio_compos_vs_seqid);
	  next;
	}
	 }

	 #####################################################################
	 ###               Getting small windows                            ##
	 #####################################################################
	 $offset = $w - int($window_size/2); # $offset starts from -5 when window_size is 10.
	 $offset = 0 if ($offset < 0);
	 if($variable_win_size ne 'v'){ $window_size = $ori_win_size;}
	 $window_1=substr($input[0], $offset, $window_size);  # window_1 is one segment
	 $window_2=substr($input[1], $offset, $window_size);  # of defined length(size)
	 @array_of_2_seq=($window_1, $window_2);

	 #####################################################################
	 ##      This is to remove the common gaps in the two windows       ##
	 #####################################################################
	 ($win_no_gap1, $win_no_gap2) = @{&main::remov_com_column(\@array_of_2_seq)};

	 #####################################################################
	 ##      Go back if the window size is too small due to gaps        ##
	 #####################################################################
	 if( (length($win_no_gap1) < $ori_win_size)&&($variable_win_size ==1) )
	 {
		$window_size+=1; $w--; next;
	 }

	 #####################################################################
	 ##      Getting Compos and Seq ids                                 ##
	 #####################################################################
	 $compos_id=${&main::compos_id_percent_array(\@array_of_2_seq)};
	 $seq_id   =${&main::seq_id_percent_array(\@array_of_2_seq)};


	 #####################################################################
	 ####    Go back if the Seq id is bigger than Compos id        #######
	 #####################################################################
	 if(($variable_win_size eq 'v') && ($seq_id >  $compos_id))
	 {
		$window_size+=1;  $w--;   next;
	 }

	 #####################################################################
	 ####   Special ID value handling                              #######
	 #####################################################################
	 if   (($seq_id    == 0  ) || ($compos_id == 0)){ $compos_id = 1;}


	 #####################################################################
	 ####     The actual calculation                               #######
	 #####################################################################
	 if(  $minus_whole_cs  eq 'm' ){  ### this substracts rates with whole CS rate
	   $ratio_compos_vs_seqid = int( $seq_id/$compos_id*10 - $ratio_whole_seq );
	   if( $ratio_compos_vs_seqid <= 0){ $ratio_compos_vs_seqid =0; }
	   if( $ratio_compos_vs_seqid > 9){ $ratio_compos_vs_seqid = 9; }
	   if( ($apply_factor eq 'f')&&($variable_win_size eq 'v') ){
		  $ratio_compos_vs_seqid =int($ratio_compos_vs_seqid * $factor); }   }
	 else{
		$ratio_compos_vs_seqid = int($seq_id/$compos_id*10);
		if( ($apply_factor eq 'f')&&($variable_win_size eq 'v') ){
		   $ratio_compos_vs_seqid =int($ratio_compos_vs_seqid * $factor); }
		if( $ratio_compos_vs_seqid > 9){ $ratio_compos_vs_seqid = 9; }
	 }

	 #$ratio_compos_vs_seqid =int(($seq_id/$compos_id)*10-$factor);}
	 #$seq_id/abs($seq_id-$compos_id+0.1)*


	 #####################################################################
	 #######       OUT of the loop (at the near to the end  ##############
	 #####################################################################
	 if( ($w + $window_size/3) > $length ){ $ratio_compos_vs_seqid='.'; }


	 #####################################################################
	 ####    When 's'how option is set(defined at prompt)          #######
	 #####################################################################
	 if( $show_calculation eq 's' ){
		printf ("SC=%-4s %-45s Seq=%-3.2s Compos=%-3.2s W=%-2s F=%-2s\n",
			  $ratio_compos_vs_seqid, $win_no_gap1,$seq_id, $compos_id, $window_size, $factor);
		printf ("        %-45s \n\n", $win_no_gap2);
	 }

	 #####################################################################
	 # Reducing increased window size according to SC rate (option 'r') ##
	 #####################################################################
	 if( ($variable_win_size eq 'v')&&($redu_window eq 'r') ){
		if( $window_size > $ori_win_size ){
	   if( $window_size > ($length/2)){ print chr(7);
	      print "\n The increased window size is over half of seq. suspicious !! \n";
	      print "\n Disable 'v' (for variable window size), at prompt and run again\n\n";
	   }
		   $window_size -= ($winsize_reduc_factor);
		   if   ($ratio_compos_vs_seqid > 7) { $winsize_reduc_factor = 3; }
		   elsif($ratio_compos_vs_seqid >  5){ $winsize_reduc_factor = 2; }
		   elsif($ratio_compos_vs_seqid >  3){ $winsize_reduc_factor = 1; }
		   elsif($ratio_compos_vs_seqid == 3){ $winsize_reduc_factor = -0.2; }
		   elsif($ratio_compos_vs_seqid == 2){ $winsize_reduc_factor = -0.4; } # This will increase the winsize
		   elsif($ratio_compos_vs_seqid == 1){ $winsize_reduc_factor = -0.8; }
		   elsif($ratio_compos_vs_seqid == 0){ $winsize_reduc_factor = -1.6; }
	}
	 }
	 push(@ratio_array, $ratio_compos_vs_seqid);
  }
  #############################################################################
  #######     FINAL outputs, 2 types                                ###########
  #############################################################################
  $ratio_whole_seq=int($ratio_whole_seq);
  return( \@ratio_array, \$ratio_whole_seq);  package main;
}


#________________________________________________________________________
# Title     : get_segment_shift_rate
# Usage     : &get_segment_shift_rate(\%hash_for_errors, \%hash_for_sec_str);
# Function  : calculates the secondary structure segment shift rate.
# Example   : <input example> First block is for the first hash input
#                             and Second is for the second hash input.
#
#             1cdg_6taa      00000442222222222242222222222777700000007000000000
#             1cdg_2aaa      00000442222222222242222222222777700000007000000000
#             2aaa_6taa      00000000000000000000000000000000000000000000000000
#
#             1cdg_6taa      -------EEE-----------EE--EEEE------EE---------EEE-
#             1cdg_2aaa      -------EEE-----------EE--EEEE------EE---------EEE-
#             2aaa_6taa      -------EEEEE------EE-EEEEEEEE----EEEE-------EEEEE-
#
#             <intermediate output example>
#             2aaa_6taa      -------00000---------00000000----0000-------00000-
#             1cdg_6taa      -------442---------------2222-----------------000-
#             1cdg_2aaa      -------222---------------2222-----------------000-
#
#             <Final output>
#             2aaa_6taa      0%
#             1cdg_6taa      67%
#             1cdg_2aaa      67%
#
# Warning   :
# Keywords  :
# Options   : 'p' or 'P' for percentage term(default)
#             'r' or 'R' for ratio term (0.0 - 1.0), where 1 means all the
#              segments were wrongly aligned.
#             's' or 'S' for Shift rate (it actually caculates the position shift
#              rate for the secondary structure segment.
#             'h' or 'H' for position Shift rate (it actually caculates the position
#              shift rate for helical segments). If this is the only option, it
#              will show the default percentage term rate for helical segments.
#              If used with 'r', it will give you ratio (0.0 - 1.0) for helical
#              segment. If used with 's' option, it will give you position shift
#              rate for only helical segments.
#             'e' or 'E' for position Shift rate (it actually caculates the position
#              shift rate for beta segments). If this is the only option, it will
#              show the default percentage term rate for beta segments. If used
#              with 'r', it will give you ratio (0.0 - 1.0) for beta. If used
#              with 's' option, it will give you position shift rate for only
#              beta segments.
# Returns   :
# Argument  : Two references of hashes. One for error rate the other for sec.
#             assignment.
# Version   : 1.1
#--------------------------------------------------------------------
sub get_segment_shift_rate{
  my($i, $k, $j, @hash, $option_string, %h, %superposed_hash,
	  $name, %out, $gap_chr, @str1, @str2, %temp, %hash_error, %hash_secondary);
  #"""""""""""""""""""""""""""""""""""""""""
  #       general argument handling        #
  #"""""""""""""""""""""""""""""""""""""""""
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^(\w)$/) ){
		  $option_string  .= $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(\w)$/) ){
		  $option_string  .= $1;    }
	  elsif(ref($_[$k]) eq "HASH") {
		  %temp = %{$_[$k]};
		  my(@keys)= sort keys (%temp);
		  my($temp_seq) = $temp{$keys[0]};

		  if($temp_seq=~/\d\d+/){
			  %hash_error = %temp; }
		  else{ %hash_secondary = %temp; }
	  }
  }#### OUTPUT are  : %hash_error  &  %hash_secondary
  #"""""""""""""""""""""""""""""""""""""""""
  #       general argument handling end    #
  #"""""""""""""""""""""""""""""""""""""""""
  %hash_secondary =%{&tidy_secondary_structure_segments(\%hash_secondary)};
  %superposed_hash =%{&superpose_seq_hash(\%hash_error, \%hash_secondary)};
  %h=%{&get_wrong_segment_rate(\%superposed_hash)};
  return(\%h);
}


#________________________________________________________________________
# Title     : search_files_in_subdir
# Usage     :
#                     $inputdir='/nfs/ind4/ccpe1/people/A Biomatic /jpo/align';
# Function  : open dir and process all files in the dir if you wish,
#             and then go in any other sub
#             if any file(dir) is linked, it skips that file.
# Example   :
# Warning   : the final var $found_from_search_files_in_subdir mustn't be 'my'ed.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : gets a ref. of a scaler (dir name) and returns nothing(void).
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub search_files_in_subdir{ package search_files_in_subdir;
   $original_dir=${$_[0]}; my($target_file)=${$_[1]};
   my(@read_files)=@{&main::read_any_dir(\$original_dir)};
   print "\n Searching for ${$_[0]}... , wait or kill me !\n";
   for $file(@read_files){ $realfile1= "$original_dir\/$file";
	 if (-l $realfile1){ next; }
	 elsif (-d $realfile1){ &main::search_files_in_subdir(\$realfile1, \$target_file); }
	 elsif (-f $realfile1){ @split =split(/\//, $realfile1); my($f) = $split[$#split];
	if($target_file eq $f){ $found_from_search_files_in_subdir =$realfile1;
	print chr(007); last;}}
	 else{ next; }  }
   return(\$found_from_search_files_in_subdir);
   last;
   package main;
}
#________________________________________________________________________
# Title     : read_any_dir
# Usage     : @file_list = @{&read_any_dir(\$absolute_path_dir_name)};
# Function  : read any dir and REMOVES the '.' and '..' entries. And then put in array.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one scaler reference.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_any_dir{
	my($in_dir, @possible_dirs, @read_files);
	if( ($_[0] eq '.') || !(defined($_[0]))){
	  $in_dir='.';
	}elsif(!(ref($_[0]))){
	  $in_dir=$_[0];
	}elsif(ref($_[0])){
	  $in_dir =${$_[0]};
	}
	if($in_dir =~ /^([\w\-\.]+)$/){
	   $in_dir="\.\/$in_dir";
	   unless(-d $in_dir){
		 $in_dir=${&dir_search_special(\$in_dir)};
	   }
	}elsif($in_dir =~/\/([\w\-\.]+)$/){
	   $in_dir="\.\/$1";  # adjust to pwd.
	   unless(-d $in_dir){
		 $in_dir=${&dir_search_special(\$in_dir)};
	   }
	}
	sub dir_search_special{
	  my($in_dir)=${$_[0]};
	  my(@ENV_dir, @probable_dir_list, @dirs,@possible_dirs, $final_dir);
	  if($in_dir =~ /\/([\w\.\-]+)$/){
		$in_dir = $1; }
	  @probable_dir_list=('ALIGN', 'PDB', 'PATH', 'HOME', 'JPO', 'PIRDIR', 'PDBSST','PDBENT',
						  'BLASTDB', 'PIRDIR', 'SWDIR');
	  for (@probable_dir_list){ @dirs=split(':', $ENV{$_});
	for (@dirs){ if (/$in_dir$/){ $final_dir = $_; } }
	  }
	  if(@possible_dirs <1){  # goes up one level and tries to find dir.
		my($pwd)=`pwd`; chomp($pwd);
	my(@temp)=split('/', $pwd);
	pop(@temp);
	my($up_pwd)=join('/', @temp);
	$in_dir="$up_pwd\/$in_dir";
	$final_dir=$in_dir if (-d $in_dir);
	  }
	  \$final_dir
	}
	opendir(DIR1,"$in_dir");
	@read_files = readdir(DIR1);
	splice( @read_files, 0, 2 );
	\@read_files;
}
#________________________________________________________________________
# Title     : scan_win_get_average (gets averages of windows of sequences of num)
# Usage     : %out1 = %{&scan_win_get_av(\%input, \$window_size, \%input2,,,,)};
#             The order of the arguments doesn't matter.
# Function  :
# Example   : input hash: ( seq1,  '13241234141234234',      (2 or more sequences accepted)
#                         seq2,  '1341324123413241234')
#             input winsize : 5;
#
#             output hash; (seq1, 1234123413241234);
#             output hash; (seq2, 1344234123412341);
#                  The numbers are ratios(compos/seqid) with given window size.
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub scan_win_get_average{
  my($window_size)=6; my($k,$r, @in_hash);
  for($k=0; $k < @_; $k++){
	if(ref($_[$k]) eq 'HASH'){ push(@in_hash, $_[$k]);  }
	elsif( ref($_[$k]) eq 'SCALAR'){ $window_size = ${$_[$k]}; }
	elsif( !(ref($_[$k])) && ( $_[$k] =~ /\d+/)){ $window_size=$_[$k]; }
  }
  for($r=0; $r < @in_hash; $r++){ my($i,$window_1, $sepa,%out_hash,$offset, $sum, $w,@win_elem_array);
	my(%input)=%{$in_hash[$r]};  my(@keys)= keys %input;
	for ($i=0; $i< @keys; $i++){ my($av);
	   $input{$keys[$i]}=~ s/,//g; $input{$keys[$i]}=~ s/[\.\- ]/0/g;
	   for($w=0; $w < length($input{$keys[$i]}); $w++){
	 $offset = $w - int($window_size/2);  if($offset < 0){ $offset = 0; }
	 $window_1= substr($input{$keys[$i]}, $offset, $window_size);
	 @win_elem_array= split(//,$window_1);
		 for(@win_elem_array){    if(/^\d[\.\d]*/){   $sum+=$_;  }    }
		 $av.=int($sum/@win_elem_array);         $sum=0;
	   }
	   $out_hash{$keys[$i]}=$av;
	}
	push(@final_out_ref, \%out_hash);
  }
  if( @final_out_ref == 1){  return($final_out_ref[0]); }
  elsif(  @final_out_ref > 1){  return(@final_out_ref); }
}

#________________________________________________________________________
# Title     : tally_2_hashes (used for get_cs_rate_for_pairs_stat.pl )
# Usage     : ($ref1, $ref2) = &tally_2_hashes(\%hash1, \%hash2, ['n', 'a', 'p', 'i']);
#              %tally_addedup=%{$ref1};    '0' position had addedup value of 1000
#              %tally_occurances=%{$ref2}; '0' position had occurred 100 times,
#                                          '0' on average had 10 in its
#                                              corresponding hash positions
# Function  : Makes hashes of tallied occurances and summed up values for disits in
#             positions.
#             calculates the occurances or occurance rates of CS rate positions.
#             The hashes should have numbers.
# Example   : you put two hash refs. (ass. array) as args (\%hash1, \%hash2)
#             The hashes are like; hash1  (name1, 0000011111, name2, 0000122222 );
#                                  hash2  (name3, 1324..1341, name4, 13424444.. );
#
#             1) The resulting 1st hash output is (0, 20,   1, 13,     2, 12)
#             which means that 0 added up to 24 in the second arg hash positions
#                              1 added up to 15 in the second arg hash positions
#                              2 added up to 18 in the second arg hash positions
#             'p' option only works with 'n' or 'a'
#             2) The resulting 2nd hash output is (0, 5,   1, 5)
#             which means that 0 occurred 5 times in the first input hash
#                              1 occurred 5 times in the first input hash
#             'p' option only works with 'n' or 'a'
# Warning   :
# Class     :
# Keywords  :  tally two hashes of numbers.
# Options   : [a n i p]
# Package   : Bio
# Reference :
# Returns   : ($ref1, $ref2), ie, two references of hash
#             averaging option causes division of 20(added up value)
#                                                by 9(occurance) in the above
#             for '0' of the first hash, so (0, 2.222,  1, 2.1666,  2, 2.4 )
#             Average is the average of numbers
#             average value in 0-9 scale (or 0-100 with 'p' option)
#             So, if there are
#                  seq1 00111110000,   The 'a' value of 0 and 1 as in the seq2
#                  seq2 33000040000    is 0-> 6/6, 1-> 4/5, while the 'n'
#                                        calc would be, 0-> 6 (60%), 1-> 4(40%)
#
# Tips      :
# Argument  : (\%hash1, \%hash2) or optionally (\%hash1, \%hash2, ['n', 'i', 'p', 'a'])
#             'n' => normalizing, 'p' => percentage out, 'i' => make int out, 'a'=> averaged
# Todo      :
# Author    : A Biomatic
# Version   : 1.2
# Used in   : get_position_shift_rate
# Enclosed  :
#--------------------------------------------------------------------
sub tally_2_hashes{
  my($factor)=10;
  my($i, $j, $k,$t,  @keys1, @keys2, %hash0, %hash1, %tally, %tally_occur,
	 %tally_all_occur, $gap_char1, $gap_char2, @string1, @string2);

   ##########################################################################################
	my(@A ) = &handle_arguments( @_ ); my( $num_opt )=${$A[7]}; my( $char_opt )=${$A[8]};
	my(@hash)  =@{$A[0]}; my(@file)   =@{$A[4]}; my(@dir   )  =@{$A[3]}; my(@array)=@{$A[1]};
	my(@string)=@{$A[2]}; my(@num_opt)=@{$A[5]}; my(@char_opt)=@{$A[6]};
   ##########################################################################################
   %hash0 = %{$hash[0]};
   %hash1 = %{$hash[1]};
   @keys1=  keys %hash0;  ### No need to sort here as you will return hash at the end
   @keys2=  keys %hash1;

  if($char_opt =~ /p/i ){ $factor =100; }

  for($i=0; $i < @keys1; $i++ ){

	 ###################################################
	 ##  Gap char detection
	 ###################################################
	 if($hash0{$keys1[$i]} =~ /([\,\-])\S+[\,\-]/){ $gap_char1 = $1; }else{ $gap_char1=''; }
	 if($hash1{$keys2[$i]} =~ /([\,\-])\S+[\,\-]/){ $gap_char2 = $1; }else{ $gap_char2=''; }


	 ###################################################
	 ##  Split the value string by gap char
	 ###################################################
	 @string1=split(/$gap_char1/, $hash0{$keys1[$i]});
	 @string2=split(/$gap_char2/, $hash1{$keys2[$i]});
	 ### @string1 => (0,0,0,0,1,1,1,1,1) @string2 => (3,4,2,13,2,1,23,3)


	 ################################################################
	 ##  Main calc part, you get %tally_all_occur and %tally_occur
	 ################################################################
	 for($j=0; $j < @string1; $j++){
	$tally_all_occur{$string1[$j]}++ ; ## <-- number of all the positions
	if( ($string2[$j]=~/[\d\^]+/)&&($string1[$j]=~/[\d\^]+/) ){
	   $tally_occur{$string1[$j]}+=$string2[$j] ; # %tally_occur is for added up counts
	}                                             # %tally_all_occur is for only the position
	 }                                                #  occurances of '0', '1' or whatever. To know
						      #  how many '0' entry were you should use this.
	 ####################################################################################
	 ##  When options were put, do more calc on %tally_all_occur and %tally_occur
	 ####################################################################################
	 if($char_opt =~ /a/i ){
	   print "\n           $char_opt ";
	   my(@cs_rates) = sort keys %tally_all_occur;
	   for($k=0; $k < @cs_rates; $k++){
		  if($tally_all_occur{$cs_rates[$k]} == 0){
			 $tally{$cs_rates[$k]} =0; next;}
		  if($char_opt =~ /i/i ){
			 $tally{$cs_rates[$k]}=int($tally_occur{$cs_rates[$k]}/$tally_all_occur{$cs_rates[$k]}); }
	  elsif($char_opt !~ /i/i ){
			 $tally{$cs_rates[$k]}= $tally_occur{$cs_rates[$k]}/$tally_all_occur{$cs_rates[$k]};
		  }
	   }
	 }
	 elsif($char_opt =~ /[np]/i){
	   my($big_sum, @cs_digits);
	   @cs_digits = sort keys %tally_occur;  # @cs_digits are (0, 1, and 2 )
	   for(@cs_digits){ $big_sum+=$tally_occur{$_}/$tally_all_occur{$_};   }
	   for($t=0; $t < @cs_digits; $t++){
		 if($big_sum ==0){ $tally{$cs_digits[$t]}=0; next; }
	 else{
	   if($char_opt =~ /i/i){
	     $tally{$cs_digits[$t]}= int(($tally_occur{$cs_digits[$t]}/$tally_all_occur{$cs_digits[$t]}/$big_sum*$factor)+0.4999);}
		   elsif($char_opt !~ /i/i ){
			 $tally{$cs_digits[$t]}= $tally_occur{$cs_digits[$t]}/$tally_all_occur{$cs_digits[$t]}/$big_sum*$factor;}
	 }
	   }
	 }
  }
  if($char_opt =~ /[an]/i){
	   print "\n           $char_opt ";
	 return(\%tally, \%tally_all_occur);}
  else{ return(\%tally_occur, \%tally_all_occur);}
}

#________________________________________________________________________
# Title     : normalize_numbers ( from 0 to 9 )
# Usage     : %output=%{&normalize_numbers(\%hash1)};
#             originally made to normalize the result of get_posi_rates_hash_out
#             in   'scan_compos_and_seqid.pl'
# Function  : with given numbers in hashes, it makes a scale of 0-9 and puts
#             all the elements in the scale. Also returns the average of the numbs.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : (\%norm_hash1, \%norm_hash2, \%norm_hash3,.... )
#
#             Example> intputhash>                   Outputhash>
#             ( '1-2', '12,.,1,2,3,4',     ( '1-2',   '9,.,0,1,2,3',
#              '2-3', '12,.,1,5,3,4',       '2-3',   '9,.,0,4,2,3',
#              '4-3', '12,3,1,2,3,4',       '3-1',   '9,3,.,.,2,3',
#              '3-1', '12,4,.,.,3,4' );     '4-3',   '9,2,0,1,2,3' );
#
# Tips      :
# Argument  : (\%hash1, %hash2, ....)
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub normalize_numbers{ my(@in)=@_;  my($split_char)=','; #<-- default split char
  my(@out_ref_of_hash, $u,$min, $max, $sum, $av, $range, @num_array,%in);
  ($min, $max, $sum, $av)=&hash_stat_for_all(@in);
  if(($max-$min)==0){ $range = 1} else { $range = ($max -$min) };
  for ($u=0; $u< @in ; $u++){ %in=%{$_[$u]};  my(@keys) = keys %in;
	 if($in{$keys[0]}=~/\,/){ $split_char=','; }  else { $split_char=''; };
	 for $name (@keys){  @num_array = split(/$split_char/, $in{$name});
		for (@num_array){   $_ = int(($_ / $range)*8) if ($_ =~ /[\-]*\d+/); }
		$in{$name}=join("$split_char", @num_array); }
	 push(@out_ref_of_hash, \%in);  }
  if( @out_ref_of_hash == 1)  {  return( $out_ref_of_hash[0]); }
  elsif( @out_ref_of_hash > 1){  return( @out_ref_of_hash   ); }
}

#________________________________________________________________________
# Title     : show_hash
# Usage     : &show_hash(\%input_hash1, \%input_hash2,.....);
# Function  : for debugging purpose. Shows any hash elems line by line in 2 columns.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub show_hash{  for($i=0; $i<= $#_; $i++){ my(%in)=%{$_[$i]};
   for(keys %{$_[$i]}){ print $_,  "\t",  $in{$_},  "\n"; }  }
}
#________________________________________________________________________
# Title     : print_seq_in_columns
# Usage     : &print_seq_in_block (\$block_leng, 'i',\%h1, 'sort', \%h2, \%hash3,,,);
# Function  : gets many refs  for one scalar  or hashes and prints
#               the contents in lines of \$block_leng(the only scalar ref. given) char.
# Example   :  With command 'print_seq_in_columns.pl c2 s2', you get:
#
#		    name1 11111111  name1 22222
#		    name2 11        name2 2222222
#		    name3 1111111   name3 22222
#		    name4 11111     name4 2222
#		    name5 11111     name5 222
#
#		    name1 3333      name1 4444
#		    name2 3333      name2 444
#		    name3 333       name3 4
#		    name4 333       name4 4444
#		    name5 3333      name5 4444444
#
# Warning   :
# Class     :
# Keywords  :
# Options   : c, i, s
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : many refs  for hash (one for bottm, one for top, etc,top hash is usually
#               to denote certain caculations or results of the bottom one
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub print_seq_in_columns{
   ###########################################################################
   my($column_num, $bl_passed, @in_ar, $c, $d,$e,$f,$g,$h,$i,$j,$k,$l, $s,$t, $x,$z,
	 $packed, $diff,$offset, $dir, $file, $in_dir, $end_found, $entry, $entry_match,
	 $gap_chr, $intl, $line, $na, $name, $larg,$names,$seq, $n_space, $offset,
	 $ordered, $output, $out_string, $pre, $pwd, $sort, $string, $tmp, $temp,
	 $title_found, $win_size, @arg_output, @in, @string_input, @k, @keys, @names,
	 @out_hash, @out_hash_final, @output_box, @outref, @read_files, @temp,
	 %hash, %input, %out_hash_final, $single_column, $largest_elem_count   );
   ##########################################################################################
	my(@A ) = &handle_arguments( @_ ); my( $num_opt )=${$A[7]}; my( $char_opt )=${$A[8]};
	my(@hash)  =@{$A[0]}; my(@file)   =@{$A[4]}; my(@dir   )  =@{$A[3]}; my(@array)=@{$A[1]};
	my(@string)=@{$A[2]}; my(@num_opt)=@{$A[5]}; my(@char_opt)=@{$A[6]};
   ##########################################################################################

  $gap_char='-';

   sub numerically{ $a <=> $b; }

  ######################################
  ####### Sub argument handling ########
  ######################################
  if($char_opt =~/n/i){ $n_space = 1; }
  if($char_opt =~/i/i){ $intl    = 1; }
  if($char_opt =~/s(\d+)/i){ $sp  = $1  }
  if($char_opt =~/c(\d+)/i){ $column_num  = $1  }
  if($char_opt =~/c[^\d]/i){ $column_num  = 1 ; } ## If just -c is given, you print single column
  if($num_opt  =~/(\d+)/i){ $bl  = $1; $bl_passed=1 if $1>0;}
  $column_num = @hash unless( defined ($column_num)  );
  $column_num = @hash if( $column_num > @hash);
  $column_num = 1 if $column_num < 1;  ## To make  '-c0' as -c1.
  if(defined($sp)){  $sp =" "x$sp }else{ $sp='  '; }


  #########  HASH input handling ############
  ## This part assigns hash keys into @names1, @names2, etc, according to the no. hash input.
  for($k=0; $k<@hash; $k++){
	   %hash=%{$hash[$k]};
	   if($sort == 1){      ## <<< this does automatic numerical handling >>>>
		  $keys_long=join("", keys %hash );
	  if ($keys_long =~ /[\d\.]+/){
	     @{"names$k"}= sort numerically keys %hash; }   # numerical sort
	  elsif($keys_long =~ /[\w\.\,]+/){
	     @{"names$k"}= sort keys(%hash);  # normal sort
		  }
	   }
	   else{
		   @{"names$k"}= keys(%hash);
		   $largest_elem_count=@{"names$k"} if (@{"names$k"} > $largest_elem_count);
	   }
	   #### IF sequence is like   'A,B,C,D,E,..' delimited by ',',  remove ','
	   for($i=0; $i< @{"names$k"}; $i++){
			  $string = $hash{${"names$k"}[$i]};
			  $larg=length($string) if $larg < length($string);
 	      if($string =~ /\,/){
		 $string=~ s/\,//g;
	      }
	  }

  }
  ## The output of above part is nothing more than  @names1, @names2, ...


  ######################################
  ####### Sub argument handling ########
  ######################################

  for($z=0; $c < @hash ; $c++){
	%hash = %{$hash[$c]};
	for($t=0;$t< @{"names$c"};$t++){
	   $s=$hash{${"names$c"}[$t]};
	   $bl = $larg if($bl_passed != 1);
	   $n = length($na)+1 if length($na) > $n;
	   if($s =~ /\-/){ $gap_char='-'; }elsif( $s =~ /\./){  $gap_char='.';  }
	   if (length($s)<$larg){
		  $offset=length($s);
		  $diff=$larg-$offset;
		  substr($s,$offset,$larg)="$gap_char"x$diff;
	   }
	}
  }
  if($intl==1){
	 $bl=$larg if ($bl_passed != 1);
	 for($k=0; $k < $larg; $k+=$bl){
		  for($c=0; $c < @hash; $c++){  # $n is the name space size
			%hash = %{$hash[$c]};
	    for($i=0; $i < @{"names$c"}; $i++){
	        $names=${"names$c"}[$i];
				$seq  =substr($hash{$names},$k,$bl);
	        printf ("%-${n}s %-${bl}s  ", $names, $seq);
	    }
	    print "\n" unless($n_space ==1);
	  }
	  print "\n";
	 }
  }
  #######################################################
  ###   This is the default printing part.
  #######################################################
  else{     # bl is the column width. should be at least 1
	 $bl=$larg if ($bl_passed != 1);
	 for($k=0; $k < $larg; $k+=$bl){
		# following is for various column number printing
		for($m=0; $m < @hash; $m+=$column_num){
		   for($i=0; $i < $largest_elem_count ; $i++){
	      for($c=$m; $c < $column_num+$m; $c++){
				 %hash = %{$hash[$c]};
				 @keys = keys (%hash);
				 $names= $keys[$i];
				 $seq=substr($hash{$names}, $k, $bl);
				 printf ("%-${n}s %-${bl}s${sp}", $names, $seq);
			  }
			  print "\n";
	   }
		   print "\n" if @names0 > 1;
		}
	 }
	 print "\n" unless($n_space ==1);
  }
}
#________________________________________________________________________
# Title     : get_position_shift_rate (derived from 'get_posi_shift_hash' )
# Usage     : %rate_hash = %{&get_position_shift_rate(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment. Takes two file names of seq.
#             Output >>
#             seq1_seq2  1110...222...2222
#             seq2_seq3  1111....10...1111
#             seq1_seq3  1111....0000.0000
#
# Example   : my(%error_rate)=%{&get_position_shift_rate(\%input, \%input2)};
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   : 'ss' for secondary structure regions(Helix and Beta region only
#                 calculation for error rate). There is specialized sub called
#              get_segment_shift_rate for sec. str. only handling.
#
#    $ss_opt            becomes    ss by  ss, SS, -ss, -SS     #  for secondary structure only
#    $H                 =         'H' by   -H or -h or H       # to retrieve only H segment
#    $S                 becomes   'S' by   -S or  S            # to retrieve only S segment
#    $E                 becomes   'E' by   -E or  E            # to retrieve only E segment
#    $T                 becomes   'T' by   -T or -t or T or t  # to retrieve only T segment
#    $I                 becomes   'I' by   -I or  I            # to retrieve only I segment
#    $G                 becomes   'G' by   -G or -g or G or g  # to retrieve only G segment
#    $B                 becomes   'B' by   -B or -b or B or b  # to retrieve only B segment
#    $HELP              becomes    1  by   -help   # for showing help
#    $simplify          becomes    1  by   -p or P or -P, p
#    $simplify          becomes    1  by   -simplify or simplify, Simplify SIMPLIFY
#    $comm_col          becomes   'C' by   -C or C or common
#    $LIMIT             becomes    L  by   -L, L               # to limit the error rate to 9 .
#
# Package   :
# Reference :
# Returns   : \%final_posi_diffs;
# Tips      :
# Argument  : %{&get_position_shift_rate(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
# Todo      :
# Author    : A Biomatic
# Version   : 1.5
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_position_shift_rate{
  #""""""""""""""""""""""< handle_arguments head Ver 1.1 >""""""""""""""""""""""""""""""""""""""
   my(@A ) = &handle_arguments( @_ ); my( $num_opt )=${$A[7]}; my( $char_opt )=${$A[8]};
   my(@hash)  =@{$A[0]}; my(@file)   =@{$A[4]}; my(@dir   )  =@{$A[3]}; my(@array)=@{$A[1]};
   my(@string)=@{$A[2]}; my(@num_opt)=@{$A[5]}; my(@char_opt)=@{$A[6]};
   my($i, $j, $c, $d, $e, $f, $g, $h, $k, $l, $p, $q, $r, $s, $t, $u, $v, $w, $x,$y,$z);
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  print __LINE__," \$char_opt is  \"$char_opt\" in get_position_shift_rate\n" if $debug eq 1;
  print __LINE__," \@string is  \"@string\" in get_position_shift_rate\n" if $debug eq 1;
  print __LINE__," \$LIMIT is  \"$LIMIT\" in get_position_shift_rate\n" if $debug eq 1;

  my(%arraySEQ)=%{$hash[0]};
  my(%arraySTR)=%{$hash[1]};
  my($gap_char, %final_posi_diffs, @stringSTR,@stringSEQ,@seq_positionSEQ,
	 @seq_positionSTR,$len_of_seq, @position_diffs, @position_corrected1,
	 @names, @whole_length, %array3, @keys_common, %DSSP_common, @stringDSSP_common);

  $gap_char='.';

  %arraySTR = %{&hash_common_by_keys(\%arraySTR, \%arraySEQ)};
  %arraySEQ = %{&hash_common_by_keys(\%arraySEQ, \%arraySTR)};
  %arraySEQ = %{&remov_com_column(\%arraySEQ)};
  %arraySTR = %{&remov_com_column(\%arraySTR)};

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if($debug eq 1){
	 print __LINE__,
	 " ## sorting sequence names. To make things constant. \n\n";  }
  @names= sort keys %arraySTR;
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  If common column of secondary structure representation option $comm_col is set
  #  open_dssp_files sub routine will get the common seq parts of all the sequences.
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if($comm_col =~ /C/i){
	 %DSSP_common=%{&open_dssp_files( @names, $H, $S, $E, $T, $I, $G, $B, $simplify, 'C')};
	 @keys_common= keys %DSSP_common;
	 @stringDSSP_common = split(/|\,/, $DSSP_common{$keys_common[0]});
	 if($debug2 eq 1){ print __LINE__," \$comm_col is set to: $comm_col \n";
		print __LINE__," \@stringDSSP_common is :@stringDSSP_common \n";
	 }
  }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  # Comparing two hashes
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for $name (@names){
	 #"""""""""""""""" Splitting the sequence string
	 if($arraySEQ{$name}=~/\,\S+\,/){
		@stringSEQ =split(/\,/, $arraySEQ{$name});
		@stringSTR=split(/\,/, $arraySTR{$name});  }
	 else{
		@stringSEQ =split(//, $arraySEQ{$name});
		@stringSTR=split(//, $arraySTR{$name});
	 }
	 print "\n",__LINE__, " \@stringSEQ  is  @stringSEQ \n" if $debug2 eq 1;
	 print "\n",__LINE__, " \@stringSTR  is  @stringSTR \n" if $debug2 eq 1;

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #   Contracting  the SEQ.
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 @seq_positionSEQ = @{&get_posi_sans_gaps(\$arraySEQ{$name})};
	 @seq_positionSTR = @{&get_posi_sans_gaps(\$arraySTR{$name})};

	 #"""""""""""""""" To get secondary structure only calc  """"""""""""""""""""""""""""
	 # It superposes the NON sec. region on  @seq_positionSTR to nullify positions.
	 #  get_posi_diff ignores non char positions in calc.
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 if( ($ss_opt =~ /ss$/i) && ($comm_col !~ /C/i) ){
	%DSSP=%{&open_dssp_files($name, $H, $S, $E, $T, $I, $G, $B, $simplify, $comm_col)};
	if($debug1 eq 1){
		  print "\n",__LINE__," open_dssp_files has options \$H ->$H \$S->$S \$E->$E \n";
		  print "\n",__LINE__," \$T->$T \$I->$I \$G->$B \$simplify->$simplify \$comm_col ->$comm_col\n";
	  &show_hash( \%DSSP );
		}
	if(ref(\%DSSP) eq 'HASH'){ # to check it %DSSP was valid, If not it skips overlaying
	   @stringDSSP = split(/|\,/, $DSSP{$name});
	   $size_of_stringDSSP = @stringDSSP;
		   $size_of_seq_positionSTR = @seq_positionSTR;
		   if($debug2 eq 1){
	       print "\n",__LINE__," \@stringDSSP is \n @stringDSSP\n";
	       print "\n",__LINE__," Size of \@stringDSSP      is $size_of_stringDSSP\n" ;
	       print "\n",__LINE__," Size of \@seq_positionSTR is $size_of_seq_positionSTR\n";
	       print "\n",__LINE__," \$gap_char is \"$gap_char\" \n" ;
		   }
		   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		   #   When the sec. str is not defined in DSSP, I delete the position of
		   #   @stringDSSP to gap(ie. make it blank to exclude error rate calc)
		   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   for($i=0; $i < @stringDSSP; $i++){
	      if($stringDSSP[$i] =~ /\W/){ $seq_positionSTR[$i]= $gap_char;}
	   }
	}
	 }elsif( $comm_col =~ /C/i){
		   print __LINE__, " Replacing position with \gap_char \"$gap_char\"\n" if $debug2 eq 1;
		   $ss_opt = 'ss'; # whether it was set or not, make it 'ss'
	   for($i=0; $i < @stringDSSP_common; $i++){
	      if($stringDSSP_common[$i] =~ /\W/){ $seq_positionSTR[$i]= $gap_char;}
	   }
	 }

	 if($debug2 eq 1){
		print __LINE__,
		print " \@seq_positionSTR is  @seq_positionSTR\n";
	 }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #   getting Position differences.
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 @position_diffs  = @{&get_posi_diff(\@seq_positionSEQ, \@seq_positionSTR)};

	 if($debug2 eq 1){
		print __LINE__,
		print " \@position_diffs is  @position_diffs\n";
	 }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  You can have two types of output according to which alignment you compare your
	 #   error rates. (1) Compare to @stringSEQ   (2) @stringSTR
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 @position_corrected1 = @{&put_position_back_to_str_seq(\@stringSEQ, \@position_diffs)};
	 $array3{$name}=join(",", @position_corrected1);

  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  # The final Step for error rate, $LIMIT is to confine error rate in one digit (ie, less than 10)
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  %final_posi_diffs =%{&get_residue_error_rate(\%array3, $LIMIT)};

  undef(@whole_length, $len_of_seq);
  return(\%final_posi_diffs);
}
#________________________________________________________________________
# Title     : get_residue_error_rate  (used in get_posi_rates_hash_out)
# Usage     : %position_diffs =%{&get_residue_error_rate(\@seq_position1, \@seq_position2)};
# Function  : This is the final step in error rate getting.
#             gets a ref. of a hash and calculates the absolute position diffs.
# Example   :
# Warning   : split and join char is ',';
# Class     :
# Keywords  :
# Options   : 'L' for limitting the error rate to 9 to make one digit output
#  $LIMIT becomes 'L' by L, l, -l, -L
# Package   :
# Reference :
# Returns   : one ref. for an array of differences of input arrays. array context.
#             ---Example input (a hash with sequences); The values are differences after
#                                comparion with structural and sequential alignments.
#             %diffs =('seq1', '117742433441...000',   <-- input (can be speparated by '' or ','.
#                      'seq2', '12222...99999.8888',
#                      'seq3', '66222...44444.8822',
#                      'seq4', '12262...00666.772.');
#             example output;
#             seq3_seq4       '0,1,0,0,0,.,.,.,,.,0,,0,0,,0,0,,.,0,,0,0,.'
#             seq1_seq2       '0,1,0,1,1,.,.,.,,.,2,,2,2,,2,2,,.,.,,2,2,1'
#             seq1_seq3       '0,1,0,1,1,.,.,.,,.,1,,1,1,,0,.,,.,.,,1,1,1'
#             seq1_seq4       '0,1,0,,1,1,.,.,.,,.,1,,1,1,0,.,.,,.,1,,2,2'
#             seq2_seq3       '0,1,0,,0,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,0'
#             seq2_seq4       '0,0,0,,1,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,.'
# Tips      :
# Argument  : Takes a ref. for hash which have positions of residues of sequences.
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   : get_position_shift_rate, previously get_each_posi_diff_hash
# Enclosed  :
#--------------------------------------------------------------------
sub get_residue_error_rate{
   my(%diffs)= %{$_[0]}; my(@names)= keys (%diffs);
   my($LIMIT)=${$_[1]} if ref($_[1]) eq 'SCALAR';
   my($LIMIT)= $_[1] unless ref($_[1]);
   my(%seqs_comp_in_pair, @temp, @temp2,$split_char, $i);
   for ($i=0; $i < @names; $i++){
	  if($diffs{$names[$i]}=~/\,/){ $split_char =',';}else{ $split_char = ''; }
	  (@{"string$i"}) = split(/$split_char/, $diffs{$names[$i]});   }
   for ($i=0; $i < @names; $i++){
	  for ($j=$i+1; $j < @names; $j ++){
	 for ($k=0; $k < @string0; $k++){
			if ((${"string$i"}[$k] =~ /[-\d+]/) && (${"string$j"}[$k] =~ /[-\d+]/)){
	       my($diff) = abs(${"string$i"}[$k] - ${"string$j"}[$k]);
			   if( ($LIMIT =~/L/i)&&($diff > 9) ){ push(@temp2, 9);
	       }else{ push(@temp2, $diff); }
			}else{ push(@temp2, '.'); } }

		 #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 #  Following if {} is for sorting output names to make  2aaa_6taa than 6taa_2aaa
		 #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 if($names[$i] <= $names[$j]){
			$seqs_comp_in_pair{"$names[$i]\_$names[$j]"}=join(",", @temp2); }
		 else{ $seqs_comp_in_pair{"$names[$j]\_$names[$i]"}=join(",", @temp2); }

		 @temp2=();
	  }
	}
   \%seqs_comp_in_pair;  # permutated
}

#________________________________________________________________________
# Title     : tell_seq_length
# Usage     : %hash_out = %{&tell_seq_length(\%hash_in)};
# Function  : tells the sequence sizes of given sequences
# Example   :
# Warning   :
# Class     : Utility
# Keywords  :
# Options   :
# Package   :
# Reference : http://sonja.acad.cai.cam.ac.uk/bioperl.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub tell_seq_length{
   ###########################################################################
   my($i,$j,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z,
	 $dir, $file, $in_dir, $end_found, $entry, $entry_match,
	 $gap_chr, $length, $line, $name, $name_found,
	 $output, $out_string, $pre, $pwd, $string,
	 @keys, @out_hash, @out_hash_final, @string, @temp,
	 @whole_file, %correct, %Final_out, %hash, %input, %out_hash_final   );
   ##########################################################################################
	my(@A ) = &handle_arguments( @II ); my( $num_opt )=${$A[7]}; my( $char_opt )=${$A[8]};
	my(@hash)  =@{$A[0]}; my(@file)   =@{$A[4]}; my(@dir   )  =@{$A[3]}; my(@array)=@{$A[1]};
	my(@string)=@{$A[2]}; my(@num_opt)=@{$A[5]}; my(@char_opt)=@{$A[6]};
   ##########################################################################################

  for($i=0; $i < @hash; $i++){
	%hash = %{$hash[$i]};
	@keys = keys %hash;
	for ($j=0; $j < @keys; $j ++){
	  if($hash{$keys[$j]}=~/\,\S+\,/){ @string= split(/\,/, $hash{$keys[$j]});
	  }else{ @string= split(//, $hash{$keys[$j]}); }
	  $h -> {$keys[$j]} = @string;  ## $h is the ref. of the anonymous hash
	}                               ## This is equivalent to "$h{$keys[$j]}= $length;"
	push(@out_hash , $h ) ;
  }
  if(@out_hash == 1){ $out_hash[0]; }
  elsif(@out_hash < 1){ die "\nSomething is wrong at tell_seq_length\n"; }
  elsif(@out_hash > 1){ return(@out_hash); }
}

#________________________________________________________________________
# Title     : handle_arguments
# Usage     : Just put the whole box delimited by the two '###..' lines below
#             to inside of your subroutines. It will call 'handle_arguments'
#             subroutine and parse all the given input arguments.
#             To use, claim the arguments, just use the variable in the box.
#             For example, if you had passed 2 file names for files existing
#             in your PWD(or if the string looks like this: xxxx.ext),
#             you can claim them by $file[0], $file[1] in
#             your subroutine.
#   #""""""""""""""""""""""< handle_arguments{ head Ver 1.2 >""""""""""""""""""""""""""""""""
#   my(@A ) = &handle_arguments( @_ ); my( $num_opt )=${$A[7]};my( $char_opt )=${$A[8]};
#   my(@hash)  =@{$A[0]};my(@file)   =@{$A[4]};my(@dir   )  =@{$A[3]};my(@array)=@{$A[1]};
#   my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};my(@raw_string)=@{$A[9]};
#   my($i, $j, $c, $d, $e, $f, $g, $h, $k, $l, $p, $q, $r, $s, $t, $u, $v, $w, $x,$y,$z);
#   if($debug==1){ print "   \@hash has \"@hash\"\n   \@raw_string has \"@raw_string\"
#   \@array has \"@array\"\n   \@char_opt has \"@char_opt\"\n   \@file has \"@file\"\n"; }
#   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#
# Function  : Sorts input arguments going into subroutines and returns default
#             arrays of references for various types (file, dir, hash, array,,,,)
#
# Example   : 'handle_arguments(\@array, $string, \%hash, 8, 'any_string')
#  ######################################################################(output example)###
#  ## $char_opt  as ('A,B,C')            ## $num_opt   as ('1,-2,3')
#  ## @char_opt  as (A, B, C)            ## @num_opt   as (1,-2, 3)
#  ## @file      as (\file1, \file2,...) ## @dir       as (\dir1, \dir2,...)
#  ## @array     as (\array1, \array2,,) ## @hash      as (\hash1, \hash2,,,,)
#  ## @string    as ('any_str_no_file_no_dir_not_opt')
#  #########################################################################################
# Warning   :
# Class     : Perl::Utility::Arg_handling
# Keywords  : handling arguments, parsing arguments,
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      : takes 0.02 u time with INDY
# Argument  : any type, any amount
# Todo      :
# Author    : A Biomatic
# Version   : 3.1
#             set_debug_option  is added.
# Used in   : everywhere
# Enclosed  :
#--------------------------------------------------------------------
sub handle_arguments{
	my($c, $d, $e, $f, $i, $j, $k, $l, $s, $t, $x, $y, $z, $char_opt, $dir, @hash,
		$file, $in_dir, $num_opt, @char_opt, @dir, @file, @string, @file_dir, @k,
		@num_opt, @string, @array, @temp, $temp);

  &set_debug_option;

  if( (@_ ==1)&& (ref($_[0]) eq 'ARRAY') ){ # when there is only 1 argument
	  @temp=@{$_[0]};
	  for $elem (@temp){
		  if( ref($temp[$elem]) ){  ## if any of the elem is a ref.
				@k = @{$_[0]};
				last;
		  }
	  }
	  @array=@{$_[0]};
  }elsif( (@_==1)&&( !ref($_[0]) ) ){
	  $temp=$_[0];
	  if(-f $temp){ push(@file, $temp) }
	  elsif(-d $temp){ push(@dir, $temp) }
	  else{ push(@string, $temp) }      ### like this &handle_arguments(\@input);
  }elsif(@_ >=1){ @k = @_ }

  #####______Start of  general argument handling______######
  for($k=0; $k < @k ;$k++){
	  if( !ref($k[$k]) ){
		  if($k[$k]=~ /^[\-]*([a-zA-Z]\d*)$/){  push(@char_opt, $1); $char_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\-([a-zA-Z]+)$/){          ## When multiple option is given,
			  my(@char_options) = split(/\,|/, $1);  push(@char_opt, @char_options);
			  $char_opt .= join("\,", @char_options); ## '-' should be used. eg. '-HEGI'
		  }elsif($k[$k]=~ /^([\-]*\d+)$/){  push(@num_opt, $1);  $num_opt  .= "$1\,";
		  }elsif(-f $k[$k]){                           push(@file,   $k[$k] );
		  }elsif(-d $k[$k]){                           push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /\/[\w\d\.\-]+[\/].+[\/]/){  push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^\/[\w\d\.\-]+[\/]*$/){     push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^[\/\w\d\-\.]+\.\w+$/){     push(@file,   $k[$k] );
		  }elsif($k[$k]=~/^\w+[\w\d\.\-]+$/){          push(@string, $k[$k] );      }
	  }elsif( ref($k[$k]) ){
		  if( ref($k[$k]) eq "SCALAR"){
			 if(${$k[$k]} =~ /^[\-]*([a-zA-Z]\d*)$/){ push(@char_opt, $1); $char_opt  .= "$1\,";
				}elsif(${$k[$k]}=~ /^\-([a-zA-Z]+)$/){ push(@char_opt, @char_options);
					$char_opt  .= join("\,", @char_options);  ## as an option string.
				}elsif(${$k[$k]}=~ /^([\-]*\d+)$/){ $num_op   .= "$1\,";  push(@num_opt, $1);
				}elsif(-f ${$k[$k]}){                            push(@file,   ${$k[$k]} );
				}elsif(-d ${$k[$k]}){                            push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\/[\/\w\d\.\-]+[\/]*$/){     push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /^[\/\w\d\-\.]+\.\w+$/){      push(@file,   ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\w+[\w\d\.\-]+$/){           push(@string, ${$k[$k]} );
				}else{                                           push(@raw_string, ${$k[$k]}); }
		  }elsif(ref($k[$k]) eq "ARRAY"){                      push(@array,  $k[$k] );
				my @temp_arr = @{$k[$k]};
				for ($i=0; $i<@temp_arr; $i++){
					if(-f $temp_arr[$i]){                            push(@file, $temp_arr[$i]);
					}elsif(-d $temp_arr[$i]){                        push(@dir , $temp_arr[$i]);
					}elsif($temp_arr[$i]=~/\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir, $temp_arr[$i] );
					}elsif($temp_arr[$i]=~/^\/[\/\w\d\.\-]+[\/]*$/){ push(@dir, $temp_arr[$i] );
					}elsif($temp_arr[$i]=~/^[\/\w\d\-\.]+\.\w+$/){   push(@file,$temp_arr[$i] );
					}elsif($temp_arr[$i]=~/^\w+[\w\d\.\-]+$/){       push(@string,$temp_arr[$i]);
					}else{                                           push(@raw_string, $temp_arr[$i]);}
				}
		  }elsif(ref($k[$k]) eq "HASH"){                       push(@hash,   $k[$k] ); }
	  }
  }
  @raw_string=(@raw_string, @string);
  @file = @{&remove_dup_in_arrayH(\@file)};
  #-----------------------------------------------------
	 sub remove_dup_in_arrayH{  my($i, @out_ref, %duplicate, @orig, @out_ref);
		for($i=0; $i<@_; $i++){  undef(%duplicate);
	 if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
	 @nondup = grep { ! $duplicate{$_}++ } @orig; push(@out_ref, \@nondup);  }
		if(@out_ref ==1){ return($out_ref[0]);}
		elsif(@out_ref >1){  return(@out_ref);}
	 }
  #-----------------------------------------------------
  return(\@hash, \@array, \@string, \@dir, \@file, \@num_opt,
			\@char_opt, \$num_opt, \$char_opt, \@raw_string );
}

#________________________________________________________________________
# Title     : set_debug_option
# Usage     : &set_debug_option;
# Function  : If you put '#' or  '##' at the prompt of any program which uses
#             this sub you will get verbose printouts for the program if the program
#             has a lot of comments.
# Example   : set_debug_option #    <-- at prompt.
# Warning   :
# Class     : Utility
# Keywords  :
# Options   : #   for 1st level of verbose printouts
#             ##  for even more verbose printouts
# $debug  becomes 1 by '#'
# $debug2 becomes 1 by '##'
#
# Package   :
# Reference : http://sonja.acad.cai.cam.ac.uk/bioperl.html
# Returns   :  $debug
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.7
#             generalized debug var is added for more verbose printouts.
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub set_debug_option{
  my($j, $i, $level);
  unless( defined($debug) ){
	for($j=0; $j < @ARGV; $j ++){
	   if( $ARGV[$j] =~/^(#+)$/){
		  print __LINE__," >>>>>>> Debug option is set by $1 <<<<<<<<<\n";
	  $debug=1; print chr(7);
		  print __LINE__," \$debug  is set to ", $debug, "\n";
	  splice(@ARGV,$j,1); $j-- ;
		  $level = length($1)+1;
		  for($i=0; $i < $level; $i++){
			 ${"debug$i"}=1;
			 print __LINE__," \$debug${i} is set to ", ${"debug$i"}, "\n";
		  }
	   }
	}
  }
}
#________________________________________________________________________
# Title     : parse_arguments  # this is the sub routine for parse_arguments.pl
# Usage     : &parse_arguments; or  (file1, file2)=@{&parse_arguments};
# Function  : Parse and assign any types of arguments on prompt in UNIX to
#             the various variables inside of the running program.
#             This is more visual than getopt and easier.
#             just change the option table_example below for your own variable
#             setttings. This program reads itself and parse the arguments
#             according to the setting you made in this subroutine or
#             option table in anywhere in the program.
# Example   :
# Warning   : HASH and ARRAY mustn't be like = (1, 2,3) or (1,2 ,3)
# Class     :
# Keywords  :
# Options   : '0'  to specify that there is no argument to sub, use
#              &parse_arguments(0);
#             parse_arguments itself does not have any specific option.
#             '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Package   :
# Reference :
# Returns   : Filenames in a reference of array
#             and input files in an array (file1, file2)=@{&parse_arguments};
# Tips      :
# Argument  : uses @ARGV
# Todo      :
# Author    : A Biomatic
# Version   : 1.2
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub parse_arguments{
  my( $c, $d, $f, $arg_num, $option_table_seen, $n, $option_filtered,
	  $option_table_example, $input_line, @input_files);

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Filtering Files from  prompt options
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( @ARGV < 1 ){ #<-- If Argument is not given at prompt
	 for(@_){
	   if($_ eq '0'){
		  last;
	   }else{
		  print "\n \"$0\" requires at least one Argument, suiciding.\n\n";
		  print chr(7); #<-- This is beeping
		  print "  To get help type \"$0  h\"\n\n\n ";
		  exit;
	   }
	 }
  }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  If there is only one prompt arg. and is [-]*[hH][elp]*, it calls
  #   &default_help and exits
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( ( @ARGV == 1 ) && ($ARGV[0] =~ /^[\-]*[hH\?][elp ]*$/) ){
	  &default_help;   exit;
  }

  for($f=0; $f < @ARGV; $f++){
	 if( ($ARGV[$f] =~ /^\w+[\-\.\w]+$/)&&(-f $ARGV[$f]) ){
	   push(@input_files, $ARGV[$f] ); next;  }
  }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   The major assignment subroutine call
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  &assign_options_to_variables;
  if($HELP == 1){ &default_help }
  return(\@input_files);
}

#________________________________________________________________________
# Title     : assign_options_to_variables
# Usage     : &assign_options_to_variables(\$input_line);
# Function  : Assigns the values set in head box to the variables used in
#             the programs according to the values given at prompt.
#             This produces global values.
# Example   : When you want to set 'a' char to a variable called '$dummy' in
#             the program, you put a head box commented line
#             '#  $dummy    becomes  a  by  -a '
#             Then, the parse_arguments and this sub routine will read the head
#             box and assigns 'a' to $dummy IF you put an argument of '-a' in
#             the prompt.
# Warning   :
# Class     :
# Keywords  :
# Options   : '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Package   : Bio::Utils
# Reference :
# Returns   : Some globaly used variables according to prompt options.
# Tips      : Used with 'parse_arguments'
# Argument  : None.
# Todo      :
# Author    : A Biomatic
# Version   : 1.4
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub assign_options_to_variables{
  my($i, $j, $op, $z, $n,);
  my($var, %val, @val, $option_table_example, @input_options);

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #      Defining small variables for option table reading
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  my($e)='=>';                  my($b)='becomes';
  my($g)='gets';                my($if)='if';
  my($is)='is';                 my($set)='is set';
  my($BY)="$by$w$wh$if";        my(@input_files);
  my($by)='by';                 my($w)='with';
  my($wh)='when';               my($o)='or';
  my(@arguments) = sort @ARGV;

  #""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Some DEFAULT $debug variables for debugging purposes
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  &set_debug_option;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   The main processing of self
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  open(SELF, "$0");    ## opens the program you ran to get the options table.
  while(<SELF>){
	 if( $first_border_and_title > 4 ){  ## This is to make it read only the first headbox.
		last;                            #  $first_border_and_title is an incremental counter.
	 }elsif( (/^ *#[_\*]{5,}$/) || (/^ *# *[Tt][itle]*[ :]*/) ){
		$first_border_and_title++;
		print __LINE__, " # assign_options_to_variables : Title line found\n" if $debug eq 1;
	 }elsif(/^ *\# *[\$\%\@].+$/){
	$op = $&;  ## $op is for the whole input option line which has $xxxx, @xxx, %xxxx format
	$op =~ s/^( *\# *)(\W\w+.+)$/$2/;  ## This is removing '#  ' in the line.
	$op =~ s/^(\W\w+.+)(\s+\#.*)$/$1/;  ## This is removing any comments in the line.

		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	  ## matching the following line input format.
	  ## $av_sc_segment     becomes    a  by  a  # To smooth the SC rates. Gets the averages of
	  ## $ARG_REG is for arguments regular expression variable. This reg. exp. matches = 'a or A or E or e' part
	  ##  which represents alternative prompt arguments possibilities. \=$b$g$is$e$set
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	  $ARG_REG = '([\-\w\.\d]+[=\w\d\.\d]*)[\;\.]* *[or\=\,\']* +([-\w\.\d]*) *[\=or\,\']* *([-\w\.\d]*) *[\=or\,\']* *([-\w\.\d]*)';

	  if($op=~/^([\$\@\%]+)([\w\-]+) *[\=$b$g$is$e$set]+ +[\'\(\[\$\@\%]*([\-\w\.\d]+)[\)\]\']* +[$by$w$wh$if]+ +$ARG_REG/){
		      ## $sym     $var        becomes            a [$a...]                           by            a -a -A

	      my $sym = $1;  #### The symbols like ($, @, %), '$' in the above.
	      my $var = $2;  #### Actual variable name 'var' from $var, 'av_sc_segment' in the above.
	      my $val = $3;  #### The becoming value  first 'a' in the above.
	      my @arg = ( $4, $5, $6, $7);  ## The alternative prompt arguments, second 'a' in the above..

			  #""""""""""""""""""""""""""""""""""""""""""""""""""""
	      #  Going through the PROMPT args.
			  #""""""""""""""""""""""""""""""""""""""""""""""""""""
	      for($z=0; $z < @arguments; $z++){     ## $arguments[$z]  is from @ARGV
		 $arguments[$z] =~ s/-//;
		 for ($i=0; $i < @arg; $i ++ ){
		   if(  ( "$arg[$i]" eq "$arguments[$z]"  ) && ($arg[$i] !~ /\=/) && ($sym eq '$') ){
		      ${"$var"}="$val";
		      splice(@arguments, $z, 1); $z --;   ## to speed up.
		      splice(@arg, $i, 1); $i --;          ## to speed up.
					  if($debug eq 1){
						 print __LINE__," \$${var} is set to \"$1\"\n";
					  }

		   }elsif( ( $arg[$i] =~ /^(\w+)\=/ )&&
						   ( $arguments[$z] =~ /^\w+\=([\w\d\.\-]+)$/)&&
						   ( $sym eq '$') ){
		      ${"$var"}="$1";
		      splice(@arguments, $z, 1); $z --;    ## to speed up.
		      splice(@arg, $i, 1); $i --;    	    ## to speed up.
					  if($debug eq 1){
						 print __LINE__,"\$${var} is set to \"$1\"\n";
					  }
		   }
		 }
	      }
	   }
	  }
   }
}
#________________________________________________________________________
# Title     : remove_dup_in_array
# Usage     : @out2 = @{&remove_dup_in_array(\@input1, \@input2,,,,)};
#             @out1 = @{&remove_dup_in_array(\@input1 )};
# Function  : removes duplicate entries in an array.
# Example   : (1,1,1,1,3,3,3,3,4,4,4,3,3,4,4);  --> (1,3,4);
# Warning   :
# Class     :
# Keywords  : merge array elements, remove_repeting_elements, remove_array_elements
# Options   :
# Package   :
# Reference :
# Returns   : one or more references.
# Tips      :
# Argument  : one or more refs for arrays.
# Todo      :
# Author    :
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remove_dup_in_array{
  my($i, @out_ref, %duplicate, @orig, @out_ref);
  for($i=0; $i<@_; $i++){
	 undef(%duplicate);
	 if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
	 @nondup = grep { ! $duplicate{$_}++ } @orig;
	 push(@out_ref, \@nondup);  }
  if(@out_ref ==1){ return($out_ref[0]);}
  elsif(@out_ref >1){  return(@out_ref);}
}
#________________________________________________________________________
# Title     : default_help         (perl version var is $] )
# Usage     : &default_help2;  usually with 'parse_arguments' sub.
# Function  : Prints usage information and others when invoked. You need to have
#             sections like this explanation box in your perl code. When invoked,
#             default_help routine reads the running perl code (SELF READING) and
#             displays what you have typed in this box.
#             After one entry names like # Function :, the following lines without
#             entry name (like this very line) are attached to the previous entry.
#             In this example, to # Function : entry.
# Example   : &default_help2; &default_help2(\$arg_num_limit);   &default_help2( '3' );
#             1 scalar digit for the minimum number of arg (optional),
#             or its ref. If this defined, it will produce exit the program
#             telling the minimum arguments.
# Warning   : this uses format and references
# Class     :
# Keywords  :
# Options   :
# Package   : File_Util
# Reference :
# Returns   : formated information
# Tips      : This usually goes with  parse_arguments.pl (= easy_opt.pl)
# Argument  :
# Todo      :
# Author    :
# Version   : 3
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub default_help{
  my($i, $perl_dir, $arg_num_limit, $pwd, $head ,$arg_num_limit );
  my($logname)=getlogin(); my($pwd)=`pwd`; my($date)=`date`; chomp($date,$pwd);
  my($not_provided)="--- not provided ---\n";
  my($file_to_read) = $0;

  for($i=0; $i < @_; $i ++){
	 if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		$arg_num_limit = ${$_[$i]};  }
	 elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		$arg_num_limit = $_[$i];     }
  }

  %entries = %{&read_head_box_help(\$file_to_read )};
  if($option_tb_found ==1){
	@option_tb=@{&read_option_table(\$file_to_read)};
  }

  foreach $help_item (keys %entries)  ## substituing with 'Not provided' message when there is no info
  {
	 ${$help_item}= $not_provided if( (${$help_item}=~/^[\W]*$/)||( !defined(${$help_item})) );
  }
  #########################################
  #########  Writing the format <<<<<<<<<<<
  #########################################
  $~ =HEADER_HELP;
  write;   ## <<--  $~ is the selection operator
  $~ =DEFAULT_HELP;
  for(sort keys %entries){  write  }
  print chr(7);  print "_"x88,"\n\n";

  if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

  if(  $option_tb_found == 1){
	#########  Printing the OPTION table contents <<<<<<<<<<<<
	print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
	   $key_press=getc();
	print @option_tb, "\n"x2 if(@option_tb > 0);
  }
format HEADER_HELP  =
_____________________________________________________________________________
			  __  __      ______      __           _____
			 /\ \/\ \    /\  ___\    /\ \         /\  _ `\
			 \ \ \_\ \   \ \ \__/    \ \ \        \ \ \L\ \
			  \ \  _  \   \ \  _\     \ \ \        \ \ ,__/
			   \ \ \ \ \   \ \ \/___   \ \ \_____   \ \ \/
				\ \_\ \_\   \ \_____\   \ \______\   \ \_\
				 \/_/\/_/    \/_____/    \/______/    \/_/  V 3
_______________________________________________________________________________
.

format DEFAULT_HELP =
 @<<<<<<<<<: @*
 $_        $entries{$_}
.
}

#________________________________________________________________________
# Title     : read_head_box_help
# Usage     : %entries = %{&read_head_box(\$file_to_read )};
# Function  : Reads the header box(the one you see on top of sub routines of
#             Jong's programs.
#             There are two types of ending line one is Jong's #---------- ...
#             the other is Astrid's  #*************** ...
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   : File_Util
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_head_box_help{
  my($title_found, %Final_out, $variable_string, $end_found, $line,
	 $entry, $entry_match );

  open(SELF, "${$_[0]}");
  my(@whole_file)=(<SELF>);
  for(@whole_file){
	if($title_found > 2){ ### This is to stop reading the file while it has found a box
	   last;              ### already.
	}elsif( / *\# *([Tt]*itle) *\: *(.*)$/){
	   $Final_out{$1}=$2;
	   $title_found ++
	}elsif( ($end_found != 1)&&($title_found==1)&&(/^\#(         +[:]* +)(.+)$/) ){
	   $Final_out{$entry_match}.= "\n $1$2";
	  # attaches to the last @entry_list element(ref)
	}elsif( ($end_found < 1)&&($title_found==1)&&(/ *\# *(\w\w\w+) *\: *(.*)/)){
	   $entry_match=$1;
	   ${"count$1"}++;
	   if( ${"count$1"} > 1){
	  $Final_out{$1}.="\n             $2";
	   }else{
		  $Final_out{$entry_match}.= $2;
	   }

	### Following is when entry line '# $certain_var = 1 by t'
	}elsif( ($end_found != 1) && ($title_found==1) && (/^\# *([\$\@\%]+.+)/) ){
	   $line = $1;
	   if($entry_match =~ /[Oo]ption/){  ## if last entry was '# Option :', attach the variable directly.
		  $Final_out{$entry_match} .= "\n             $line";
	   }else{                            ## if last entry wasn't '# Option :', find Option
		  for $entry (keys %Final_out){  ##  and attach the variable to it
			 if ($entry =~ /[Oo]ption/){
				$Final_out{$entry} .= "\n             $line";
			 }
		  }
	   }
	}elsif( ($title_found==1)&&(/ *\#[\*\-]{12,}/)){  ## to match '#-----..' or '#*******..'(Astrid's)
	   $end_found++;
	}elsif( (/^#{10,} option table of this program   #{10,}/)&&($end_found >=1) &&($title_found==1)){
	   $option_tb_found++; ### This is a global var.
	}
  }
  print "\n $option_tb_found  E $end_found T $title_found\n";
  \%Final_out;
}                  ##<<--- ENd of the sub read_head_box

#________________________________________________________________________
# Title     : read_option_table
# Usage     :
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   : File_Util
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_option_table{
   my($table_found, @option_tb, $head);
	open(SELF, "${$_[0]}");
	while(<SELF>){
	  if( (/^ *#+/) && ( $table_found== 1) ){
		push (@option_tb, "$_");
	  }elsif( ($table_found != 1)&&(/^ *\#+ *[Oo]ption *[Tt]able */) ){
		 $table_found=1; $head="############## Option Table for $logname\'s \"$0\"\n"; ##
		 push(@option_tb, $head);
	  }
	  if( ($table_found==1)&&(/^ *###################+ *$/)){  ### to find the end point of reading
		 $table_found =0; last; }
	}
	\@option_tb;
}





#________________________________________________________________________
# Title     : hash_stat_for_all
# Usage     : %out=%{&hash_average(\%in, \%in2,..)};
# Function  : gets the min, max, av, sum for the whole values of ALL the
#             hashes put in. (grand statistics)
# Example   : %in =(1, "13242442", 2, "92479270", 3, "2472937439");
#             %in2=(1, "28472", 2, "23423240", 3, "123412342423439");
#
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : normal array of ($min, $max, $sum, $av)
#             Example out
#                -----------------------------------
#                of the whole    |   0   9  110   6
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub hash_stat_for_all{ my(@out_av_hash, $v,@num_arr,$elem,$sum,$min, $num_all,$max,$split_char);
  for($v=0; $v<@_; $v++){ my(%input)=%{$_[$v]};
	 for $name(keys(%input)){
	   if($input{$name} =~ /\,/){ $split_char=','; }
	   else{ $split_char=''; } @num_arr=split(/$split_char/, $input{$name});
	   for $elem(@num_arr){
	 if($elem =~/[\-]*\d+/){   $min=$elem unless(defined($min));
	    $min =$elem if $elem < $min; $max =$elem if $elem > $max;
	    $sum+=$elem; $num_all++; } } } }
  if($num_all == 0){ $av=0; $sum=0; $min=0; $max=0; }else { my($av)=$sum/$num_all; }
  push(@out_array, ($min, $max, $sum, $av));
  return(@out_array); ## not a ref. for an array !!!
}

#________________________________________________________________________
# Title     : get_posi_rates_hash_out_compact (derived from 'get_posi_shift_hash' )
# Usage     : %rate_hash = %{&get_posi_shift_hash(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment.
#             Output >>  something like below but, without gaps, so final one is;
#             seq1_seq2  1110...222...2222     seq1_seq2  11102222222
#             seq2_seq3  1111....10...1111  -> seq2_seq3  1111101111
#             seq1_seq3  1111....0000.0000     seq1_seq3  111100000000
#
# Example   :
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : \%final_posi_diffs_compact;  Compare with  'get_posi_rates_hash_out_jp'
# Tips      :
# Argument  : %{&get_posi_rates_hash_out(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_rates_hash_out_compact{ my(%array1)=%{$_[0]};  my(%array2)=%{$_[1]};
  my(@string1,@string2,@seq_position1,@seq_position2,$len_of_seq,@position_diffs,
	 @position_corrected1,@names, @whole_length, %array3);
  %array1 = %{&hash_common_by_keys(\%array1, \%array2)};
  %array2 = %{&hash_common_by_keys(\%array2, \%array1)};
  %array1 = %{&remov_com_column(\%array1)};
  %array2 = %{&remov_com_column(\%array2)};
  @names= keys %array2;
  for $name (@names){
	 @string1 =split('', $array1{$name});  @string2 =split('', $array2{$name});
	 @seq_position1 = @{&get_posi_sans_gaps(\$array1{$name})};
	 @seq_position2 = @{&get_posi_sans_gaps(\$array2{$name})};
	 $len_of_seq =(@seq_position2);  push(@whole_length, $len_of_seq);
	 @position_diffs = @{&get_posi_diff(\@seq_position1, \@seq_position2)};
	 $array3{$name}=join(",", @position_diffs);  }
  my(%final_posi_diffs_compact)=%{&get_each_posi_diff_hash(\%array3)}; undef(@whole_length, $len_of_seq);
  return(\%final_posi_diffs_compact);
}

#________________________________________________________________________
# Title     : get_posi_rates_hash_out_jp (derived from 'get_posi_shift_hash' )
# Usage     : %rate_hash = %{&get_posi_shift_hash(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment. get_posi_rates_hash_out_jp
#             results in jp template sequence, while get_posi_rates_hash_out_msf does
#             in msf template sequence.
#             Output >>
#             seq1_seq2  1110...222...2222   <--- the alignment template is JPO's
#             seq2_seq3  1111....10...1111        (ie structural)
#             seq1_seq3  1111....0000.0000
#
# Example   :
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : \%final_posi_diffs;
# Tips      :
# Argument  : %{&get_posi_rates_hash_out_jp(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_rates_hash_out_jp{ my(%array1)=%{$_[0]};  my(%array2)=%{$_[1]};
  my(@string1,@string2,@seq_position1,@seq_position2,$len_of_seq,@position_diffs,
	 @position_corrected1,@names, @whole_length, %array3);
  %array1 = %{&hash_common_by_keys(\%array1, \%array2)}; %array2 = %{&hash_common_by_keys(\%array2, \%array1)};
  %array1 = %{&remov_com_column(\%array1)};              %array2 = %{&remov_com_column(\%array2)};
  @names= keys %array2;
  for $name (@names){
	 @string1 =split('', $array1{$name});  @string2 =split('', $array2{$name});
	 @seq_position1 = @{&get_posi_sans_gaps(\$array1{$name})}; @seq_position2 = @{&get_posi_sans_gaps(\$array2{$name})};
	 $len_of_seq =(@seq_position2);  push(@whole_length, $len_of_seq);
	 @position_diffs = @{&get_posi_diff(\@seq_position1, \@seq_position2)};
	 @position_corrected1 = @{&put_position_back_to_str_seq(\@string2, \@position_diffs)};
	 $array3{$name}=join(",", @position_corrected1);  }
  my(%final_posi_diffs)=%{&get_each_posi_diff_hash(\%array3)}; undef(@whole_length, $len_of_seq);
  return(\%final_posi_diffs);
}



#________________________________________________________________________
# Title     : get_posi_rates_hash_out (derived from 'get_posi_shift_hash' )
# Usage     : %rate_hash = %{&get_posi_shift_hash(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment.
#             Output >>
#             seq1_seq2  1110...222...2222
#             seq2_seq3  1111....10...1111
#             seq1_seq3  1111....0000.0000
#
# Example   :
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : \%final_posi_diffs;
# Tips      :
# Argument  : %{&get_posi_rates_hash_out(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_rates_hash_out{
  my(%array1)=%{$_[0]};
  my(%array2)=%{$_[1]};
  my(@string1, @string2, @seq_position1, @seq_position2,
	 $len_of_seq,@position_diffs, @position_corrected1,
	 @names, @whole_length, %array3);
  %array1 = %{&hash_common_by_keys(\%array1, \%array2)};
  %array2 = %{&hash_common_by_keys(\%array2, \%array1)};
  %array1 = %{&remov_com_column(\%array1)};
  %array2 = %{&remov_com_column(\%array2)};
  @names= keys %array2;
  for $name (@names){
	 @string1 =split('', $array1{$name});
	 @string2 =split('', $array2{$name});
	 @seq_position1 = @{&get_posi_sans_gaps(\$array1{$name})};
	 @seq_position2 = @{&get_posi_sans_gaps(\$array2{$name})};
	 $len_of_seq =(@seq_position2);
	 push(@whole_length, $len_of_seq);
	 @position_diffs = @{&get_posi_diff(\@seq_position1, \@seq_position2)};
	 @position_corrected1 = @{&put_position_back_to_str_seq(\@string2, \@position_diffs)};
	 $array3{$name}=join(",", @position_corrected1);  }
  my(%final_posi_diffs)=%{&get_each_posi_diff_hash(\%array3)};
  undef(@whole_length, $len_of_seq);
  return(\%final_posi_diffs);
}

#________________________________________________________________________
# Title     : get_posi_rates_hash_out (derived from 'get_posi_shift_hash' )
# Usage     : %rate_hash = %{&get_posi_shift_hash(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment.
# Example   :
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : \%final_posi_diffs;
# Tips      :
# Argument  : %{&get_posi_rates_hash_out(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
#             Output >>
#             seq1_seq2  1110...222...2222
#             seq2_seq3  1111....10...1111
#             seq1_seq3  1111....0000.0000
#
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_rates_hash_out_msf{
  my(%array1)=%{$_[0]};
  my(%array2)=%{$_[1]};
  my(@string1, @string2, @seq_position1, @seq_position2,
	 $len_of_seq,@position_diffs, @position_corrected1,
	 @names, @whole_length, %array3, $name);
  %array1 = %{&hash_common_by_keys(\%array1, \%array2)};
  %array2 = %{&hash_common_by_keys(\%array2, \%array1)};
  %array1 = %{&remov_com_column(\%array1)};
  %array2 = %{&remov_com_column(\%array2)};
  @names= keys %array2;
  for $name (@names){
	 @string1 =split(/\,|/, $array1{$name});
	 @string2 =split(/\,|/, $array2{$name});
	 @seq_position1 = @{&get_posi_sans_gaps(\$array1{$name})};
	 @seq_position2 = @{&get_posi_sans_gaps(\$array2{$name})};
	 $len_of_seq =(@seq_position1);
	 push(@whole_length, $len_of_seq);
	 @position_diffs = @{&get_posi_diff(\@seq_position1, \@seq_position2)};
	 @position_corrected1 = @{&put_position_back_to_str_seq(\@string2, \@position_diffs)};
	 $array3{$name}=join(",", @position_corrected1);
  }
  my(%final_posi_diffs)=%{&get_each_posi_diff_hash(\%array3)};
  undef(@whole_length, $len_of_seq);
  # show_hash(\%final_posi_diffs);
  return(\%final_posi_diffs);
}


#________________________________________________________________________
# Title     : compos_id_percent_array  (more than 2 elements array)
# Usage     : $percent = &compos_id_percent_array(@any_array_sequences);
#             The way identity(composition) is derived is;
# Function  : produces amino acid composition identity of any given number of sequences.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   : get_windows_compos_and_seqid_rate_array,
#
# Enclosed  :
#--------------------------------------------------------------------
sub compos_id_percent_array{ my(@input)=@{$_[0]}; my($largest,$iden,@temp,$iden_sum,$final_iden, @all_pairs_id);
   for($i=0; $i<=$#input; $i++){  $input[$i]=~ tr/a-z/A-Z/;  $input[$i]=~ s/[\.\-\s]//g;
	  @temp = split('', $input[$i]);  (@{"string$i"})= @temp;
	  $largest = @{"string$i"} if @{"string$i"} > $largest;    }
   for($i=0; $i< @input; $i++){ #_________ permutating ___________
	  for($j=$i+1; $j<=$#input; $j++){
		 for ($k=0; $k <= $largest; $k ++){  # getting composition tables for two seqs.
			$compos_table1{${"string$i"}[$k]}++ if (${"string$i"}[$k] =~ /\w/);
	    $compos_table2{${"string$j"}[$k]}++ if (${"string$j"}[$k] =~ /\w/);   }
	 $iden =${&calc_compos_id_hash(\%compos_table1, \%compos_table2)};
	 push(@all_pairs_id, $iden);  %compos_table1=();  %compos_table2=();   }   }
   for $iden (@all_pairs_id){  $iden_sum+=$iden;  }
   $final_iden=$iden_sum/(@all_pairs_id);
   #-----------------------------------------------------
   #  Input here is like :  %hash1= (A,3,B,3,C,4,D,4), %hash2= (A,4,B,1,C,4)
   sub  calc_compos_id_hash{  # input is like this;
	 my(%hash1)=%{$_[0]}; my(%hash2)=%{$_[1]}; my(%common_of_the_2)=();
	 my($common, $compos_id, $sum_residues, $sum_of_the_common_residue_no);
	 my(@values1) = values (%hash1);   my(@values2) = values (%hash2);
	 my(@combined_values)=(@values1,@values2);
	 for $elem (@combined_values){  $sum_residues += $elem;   }
	 if($sum_residues == 0){ $compos_id =0; } # to prevent Illegal division error.
	 else{ for $key1(keys %hash1){
	     $common=&smaller_one($hash1{$key1}, $hash2{$key1});
	        sub smaller_one{if ($_[0] > $_[1]){ $_[1];}else{$_[0];}}
	     $sum_of_the_common_residue_no += $common;     }
		  $compos_id = $sum_of_the_common_residue_no/($sum_residues/2)*100;   }
	 \$compos_id;
   }
   #-----------------------------------------------------
   return ( \$final_iden ); # final identity for any given set of strings(seq).
}

#________________________________________________________________________
# Title     : seq_id_percent_array  (more than 2 elements array)
# Usage     : $percent = &seq_id_percent_array(@any_array_sequences);
#             The way identity(pairwise) is derived is;
#
# Function  : produces amino acid composition identity of any given number of sequences.
# Example   :
# Warning   : This can handle 'common gaps' in the sequences
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub seq_id_percent_array{ my(@input, $denominator,@all_pairs_id, $percent_id);
   my($largest,$p,$i,$j,$k,$iden_residue_num,$iden,@temp,$iden_sum,$gap_num,$final_iden);
   for($d=0; $d<@_; $d++){
	  if(ref($_[$d]) eq 'ARRAY'){ @input=@{$_[$d]}; }
	  elsif( (ref($_[$d]) eq 'SCALAR') &&( ${$_[$d]}=~/^[aA]/) ){ $average_len_opt =1 }
	  elsif( !(ref($_[$d])) && ( $_[$d] =~/^[aA]/) ){ $average_len_opt =1 } }
   if ((@input== 1)||( @input== 0)){
	  print "\n\n \" $0 \"  requires at least 2 sequences\n\n";
	  print "\n Abnormally dying at seq_id_percent_array in $0 \n\n"; print chr(7); exit;}
   $shortest=length($input[0]); my($sans_gap_seq, $length_sum, $average_seq_len);
   for($p=0; $p < @input; $p++){
	  $input[$p]=~ tr/a-z/A-Z/; $sans_gap_seq=$input[$p];  $sans_gap_seq=~s/\W//g;
	  $input[$p]=~ s/\W/./g;  (@{"string$p"})=split('', $input[$p]);
	  $largest = length($input[$p]) if length($input[$p]) > $largest;
	  $shortest = length($sans_gap_seq) if length($sans_gap_seq) < $shortest;
	  $length_sum += length($sans_gap_seq);     }
   $average_seq_len = $length_sum/@input;
   for($i=0; $i< @input; $i++){
	  for($j=$i+1; $j< @input; $j++){
	 for ($k=0; $k <  $largest; $k ++){  # getting composition tables for two seqs.
	    if ((${"string$i"}[$k] !~ /\W/)&&(${"string$i"}[$k] eq ${"string$j"}[$k])){
			   $iden_residue_num++; }
	    elsif((${"string$i"}[$k] =~ /\W/)&&(${"string$i"}[$k]=~ /\W/)){ $gap_num++; }}
	 if( $average_len_opt == 1){ $denominator = $average_seq_len; }
	 else{ $denominator = $shortest; }
	 if($denominator == 0){ $denominator=1; }  # in the above it is 50% rather than 0.07%
	 $percent_id=($iden_residue_num/($denominator))*100;
	 push(@all_pairs_id, $percent_id);
	 undef ($iden_residue_num, $gap_num); } }
   for (@all_pairs_id){  $iden_sum+=$_;    }
   $final_iden=$iden_sum/($#all_pairs_id+1); return( \$final_iden );
}

#________________________________________________________________________
# Title     : min_elem_array
# Usage     : ($out1, $out2)=@{&min_elem_array(\@array1, \@array2)};
#             ($out1)       =${&min_elem_array(\@array1)          };
# Function  : gets the smallest element of any array of numbers.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one or more ref. for scalar numbers.
# Tips      :
# Argument  : numerical arrays
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub min_elem_array{  my(@out_min_elem, @input, $min_elem);
  for($i=0; $i< @_; $i++){ @input=@{$_[$i]}; $min_elem=$input[$#input];
	for (@input){  $min_elem=$_ if ((/[\-\d]+/)&&($_ < $min_elem));    }
	push(@out_min_elem, $min_elem);}
  if(@_ == 1){ return( \$min_elem ); }
  elsif(@_ > 1 ){  return( \@out_min_elem ) };
}

#________________________________________________________________________
# Title     : open_dssp_files
# Usage     : (*out, *out2) = @{&open_dssp_files(\$inputfile1, \$inputfile2, \$H, \$S,,,,)};
#             (@out)        = @{&open_dssp_files(\$inputfile1, \$inputfile2, \$H, \$S,,,,)};
# Function  : open dssp files and put sequences in a hash(s)
#              It can take options for specific secondary structure types. For example,
#              if you put an option $H in the args of the sub with the value of 'H'
#              open_dssp_files will only read secondary structure whenever it sees 'H'
#              in xxx.dssp file ignoring any other sec. str. types.
#              If you combine the options of 'H' and 'E', you can get only Helix and long
#              beta strand sections defined as segments. This is handy to get sec. str. segments
#              from any dssp files to compare with pdb files etc.
#             With 'simplify' option, you can convert only all the 'T', 'G' and 'I' sec. to
#              'H' and 'E'.
# Example   :
# Warning   : 6taa.dssp  and 6taa are regarded as the same.
# Class     :
# Keywords  :
# Options   : H, S, E, T, I, G, B, P, C, -help
# $H        =        'H' by   -H or -h or H or h  # to retrieve 4-helix (alpha helical)
# $S        becomes  'S' by   -S or -s or S or s  # to retrieve Extended strand, participates in B-ladder
# $E        becomes  'E' by   -E or -e or E or e  # to retrieve residue in isolated Beta-bridge
# $T        becomes  'T' by   -T or -t or T or t  # to retrieve H-bonded turn
# $I        becomes  'I' by   -I or -i or I or i  # to retrieve 5-helix (Pi helical) segment output
# $G        becomes  'G' by   -G or -g or G or g  # to retrieve 3-helix (3-10 helical)
# $B        becomes  'B' by   -B or -b or B or b  # to retrieve only B segment
# $simplify becomes   1  by   -p or P or -P, p
# $comm_col becomes  'c' by   -c or c or C or -C or common
# $HELP     becomes   1  by   -help   # for showing help
#
# Package   :
# Reference :
# Returns   : (*out, *out2)  or (@out_array_of_refs)
# Tips      :
# Argument  : files names like (6taa, 6taa.dssp) If you put just '6taa' without extension, it
#             searches if there is a '6taa.dssp' in both PWD and $DSSP env. set directory.
#             ---------- Example of dssp ---
#             **** SECONDARY STRUCTURE DEFINITION BY THE PROGRAM DSSP, VERSION JUL
#             REFERENCE W
#             HEADER    RIBOSOME-INACTIVATING PROTEIN           01-JUL-94   1MRG
#             COMPND    ALPHA-MOMORCHARIN COMPLEXED WITH ADENINE
#             SOURCE    BITTER GOURD (CUCURBITACEAE MOMORDICA CHARANTIA) SEEDS
#             AUTHOR    Q
#             246  1  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
#             112 95.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .
#             171 69.5   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .
#             12   4.9   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
#             36  14.6   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
#             1    0.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
#             1    0.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
#             74  30.1   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
#             5    2.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
#             1    2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           .
#             0    0  0  0  1  1  0  2  0  0  1  0  0  1  0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0    RESIDUES PER ALPHA HELIX         .
#             1    0  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    PARALLEL BRIDGES PER LADDER      .
#             2    0  1  2  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    ANTIPARALLEL BRIDGES PER LADDER  .
#             2    0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
#             #   RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
#             1    1   D              0   0  132    0, 0.0   2,-0.3   0, 0.0  49,-0.2   0.000 360.0 360.0 360.0 153.4   44.0   96.9  -23.8
#             2    2   V  E     -a   50   0A  10   47,-1.5  49,-2.8   2, 0.0   2,-0.3  -0.889 360.0-163.3-115.9 151.4   43.1  100.4  -22.5
#             3    3   S  E     -a   51   0A  63   -2,-0.3   2,-0.3  47,-0.2  49,-0.2  -0.961  10.3-172.8-131.0 152.3   44.8  103.7  -23.4
#             4    4   F  E     -a   52   0A   8   47,-2.2  49,-2.3  -2,-0.3   2,-0.4  -0.985   6.9-161.2-143.2 139.5   45.0  107.2  -22.0
#             5    5   R  E     -a   53   0A 144   -2,-0.3   4,-0.2  47,-0.2  49,-0.2  -0.993   9.7-156.0-121.0 125.9   46.6  110.2  -23.6
#             6    6   L  S    S+     0   0    1   47,-2.3   2,-0.5  -2,-0.4   3,-0.4   0.644  73.2  90.9 -73.3 -22.4   47.5  113.2  -21.4
#             7    7   S  S    S+     0   0   81   47,-0.3   3,-0.1   1,-0.2  -2,-0.1  -0.695 106.0   5.2 -75.5 121.0   47.4  115.6  -24.4
#             8    8   G  S    S+     0   0   72   -2,-0.5  -1,-0.2   1,-0.3   5,-0.1   0.269  97.6 147.8  90.2 -10.7   43.9  117.0  -24.7
#             9    9   A        +     0   0   10   -3,-0.4  -1,-0.3  -4,-0.2  -3,-0.1  -0.256  16.8 166.8 -58.8 142.4   42.9  115.2  -21.5
#             (\$inputfile1, \$inputfile2, .... )};
# Todo      :
# Author    : A Biomatic
# Version   : 2.8
#             $debug feature has been added to make it produce error messages with '#' option.
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_dssp_files{
  #""""""""""""""""""""""< handle_arguments head Ver 1.0 >""""""""""""""""""""""""""""""""""""""
   my($i, $j, $c, $d, $e, $f, $g, $h, $k, $l, $p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z,
	 $dir, $file, $line, $name, $ori_name, @keys, @names,  $str_input, $C, @out_hash_ref_list,
	 @temp, $temp, %hash, %input, %com_col_hash, $ref_hash_out, $BASE );
   my(@A ) = &handle_arguments( @_ ); my( $num_opt )=${$A[7]}; my( $char_opt )=${$A[8]};
   my(@hash)  =@{$A[0]}; my(@file)   =@{$A[4]}; my(@dir   )  =@{$A[3]}; my(@array)=@{$A[1]};
   my(@string)=@{$A[2]}; my(@num_opt)=@{$A[5]}; my(@char_opt)=@{$A[6]};
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  if($char_opt !~ /[HEBGIST]/i){  ## This is default sec. str type setting. (full representation)
	  $char_opt = 'HEBGIST';
  }
  if ($debug eq 1){
	 print __LINE__, " # open_dssp_files : \$simplify     is  $simplify\n" ;
	 print __LINE__, " # open_dssp_files : \@file given   is  @file \n" ;
	 print __LINE__, " # open_dssp_files : \@string given is  @string\n" ;
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ### Big main loop for input argument handling   ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  This is to check if the given file is not in pwd but in ENV var $DSSP
  #  Or if the file name was given only by the base name of seq(eg. 1cdg rather
  #  than 1cdg.dssp
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for ($i=0; $i < @string; $i ++){
	   print __LINE__, " ${i}th string input is  $string[$i] \n" if $debug eq 1;
	   $string[$i] = "$string[$i]\.dssp"; ## adding  .dssp extension
	   print __LINE__, " ${i}th string inputwith \.dssp is now, $string[$i] \n" if $debug eq 1;

	   if(-f $string[$i]){
		  print chr(7) if $debug eq 1;
		  print __LINE__, " Your input filename exist in this File: $string[$i]\n" if $debug eq 1;
		  unshift(@file, "$string[$i]");
	   }
	   elsif(-l $string[$i]){
		  print chr(7) if $debug eq 1;
		  print "\n Your input filename exist as a Link to : $string[$i]\n" if $debug eq 1;
		  unshift(@file, "$string[$i]");
	   }
	   elsif( -d $ENV{'DSSP'} ){
		  $string[$i] =~ s/(\w+)\.dssp$/$1/; ## stripping .dssp extension
		  if( -e "$ENV{'DSSP'}\/$string[$i]\.dssp" ){
			unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
			$BASE = $string[$i];
		  }else{
			 print chr(7);
			 print __LINE__, " !! Error your DSSP env setting seems wrong. \n";
			 print __LINE__, " !! Your DSSP env path is also a link. \n" if (-l $ENV{'DSSP'});
			 print __LINE__, " I can't find  $ENV{'DSSP'}\/$string[$i] \n\n";
		  }
	   }
	   elsif( -l $ENV{'DSSP'} ){ #"""""""  IF $DSSP was a link
		  print __LINE__, " !! Your DSSP env path is also a link. \n" if $debug eq 1;
		  if( -e "$ENV{'DSSP'}\/$string[$i]\.dssp" ){
			unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
			$BASE = $string[$i];
		  }
		  elsif( -e "$ENV{'DSSP'}\/$string[$i]" ){
			unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
		  }
	   }
  }

  @file=@{remove_dup_in_array(\@file)};

  if ($debug eq 1){
	 print __LINE__, " # open_dssp_files : ENV set for dssp is $ENV{'DSSP'} \n" ;
	 print __LINE__, " # open_dssp_files : Final \@file given are \" @file \"\n" ;
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  END of File and string input checking in searching for the right dssp file.
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"""""""""""""""""""""""" MAIN """"""""""""""""""""""""""""""""""""""
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for($i=0; $i< @file; $i++){  ## <<-- loops over the input files.
	   my($flag, %hash, $name, $s, $matched, $ori_name, $chain);
	   my($real_file) = $file[$i];
	   $file[$i] =~ s/(.*\/)(\w+)\.(\w+)$/$2/; ## stripping .dssp extension
	   $file[$i] =~ s/(\w+)\.(\w+)$/$1/;       ## stripping .dssp extension
	   $ori_name = $name = $file[$i];
	   print "\n",__LINE__, " VAR \$ori_name is  $ori_name , \$file\[\$i\] is $file[$i]\n" if $debug eq 1;
	   unless(-e $real_file){
		 print "\n",__LINE__,"  !!! ERROR $real_file does not exists as the final filename\n" if $debug eq 1;
		 splice(@file, $i, 1); $i--;
		 print "\n",__LINE__,"  Skipping to the next file to open" if $debug eq 1;
		 next;
	   }

	   open(FILE_1,"$real_file");
	   print "\n",__LINE__, " ${i}th file $real_file is being opened from \@file \n" if $debug eq 1;
	   print "_"x86,"\n", if $debug eq 1;

	   while(<FILE_1>){
		  if(/^[\s]*\#\s+RESIDUE/){
			 $flag =1;
			 print __LINE__," \"#  RESIDUE\"   string found at line $. in $real_file\n" if $debug eq 1;
			 next
		  }  ##   '#  RESIDUE' is the starting key

		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  #    Matching the column
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	  if(  ($flag==1) && (/^[\s]*-*\d+\s+-*\d*\s+[\w]\s\s([\w ]) /)  ){
			 $matched = $1;
			 print __LINE__," \"$matched\" is matched\n" if $debug2 eq 1;

			 if( $char_opt =~ /$matched/){ ## Here OPTIONS affect the operation.
				$s .= $matched;    ## $match_option is like 'HE'. If the
				next;              ## single char $matched is H or E, it will be
			 }else{                ## annexed to $s as an output.
				$s .= ' ';
				next;
			 }  # <-- this is necessary to get the right length (not to ignore
		  }     #     not matched char by converting them to ' '.

		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  #    When there are chains like A, B, ,,,
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  elsif( ($flag==1) && (/^[\s]*-*\d+\s+-*(\d+)[\s\w]+(\w)\s+[\w]\s\s([\w ]) /) ){
			 $chain = $2;   ## $flag  is for the starting key
			 # ${"chain_start$name$2"} = $1 unless defined(${"chain_start$name$2"});
			 my($matched_chain) = $3;
			 if( $char_opt =~ /$matched_chain/){
			   $s .= $matched_chain;   next; }
			 else{
			   $s .= ' '; next; }
		  }elsif( (/^\s+\d+\s+\!/)&&($chain =~/\w/) ){
			 $name="$name$chain";
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ##   IF simplify  option is set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 if($simplify eq 1){
				$s =~ tr/TGI/EHH/;   ### change the characters.
				print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
			 }
			 if($debug eq 1){ print __LINE__, " Name of seq:  $name \n"; }
			 $hash{$name}=$s; $s='';
			 $name=$ori_name; next;
		  }
	   }
	   close(FILE_1);  ##<<---- Reading finished.

	   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   ##  Naming procedure
	   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   if($chain =~/^\w$/){     # when there are chains, put A,B, etc to seq. names.
		  $name="$name$chain";  ## <<-- This is for the last chain entry.
		  if($debug eq 1){ print __LINE__, " Name of seq:  $name, There were Chains !\n"; }
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  ##   IF simplify  option is set
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  if($simplify eq 1){
			 $s =~ tr/TGI/EHH/;   ### change the characters.
			 print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
		  }
		  $hash{$name}=$s;
		  $s='';   ##<<--- This is essential, a former bug!
		  $name=$ori_name;
	   }else{      # <<-- Without chains option.
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  ##   IF simplify  option is set
		  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		  if($simplify eq 1){
			 $s =~ tr/TGI/EHH/;   ### change the characters.
			 print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
		  }
		  $hash{$name}=$s;
		  if($debug eq 1){ print __LINE__, " Name of seq:  $name \n"; }
		  $s='';   #<<--- This is essential, a former bug!
	   }

	   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   ### OUTput format determination according to options #####
	   #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   if($debug eq 1){ print "\n", __LINE__, " The Hash out of \"$real_file\" is \n ";
		  &show_hash(%hash);
	   }
	   push(@out_hash_ref_list, \%hash) if ref(\%hash) eq 'HASH';
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"""""""""""" END of Main """""""""""""""""""""""""""
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  if($comm_col =~ /c/i){ # $comm_col  is a global
	 if($debug eq 1){
		print "\n", __LINE__;
		print " # open_dssp_files : you have put 'c' option for common column only\n";
		$temp = @out_hash_ref_list;
		print __LINE__, " # open_dssp_files : No. of hashes passed to get_common_column is: $temp\n";
		print __LINE__, " # open_dssp_files : The hash are(is) : @out_hash_ref_list\n";
	 }
	 $ref_hash_out = &get_common_column(@out_hash_ref_list);
	 return($ref_hash_out);
  }else{
	 if(@out_hash_ref_list == 1){ return($out_hash_ref_list[0]); }
	 elsif(@out_hash_ref_list > 1){ return(@out_hash_ref_list);  }
  }
}

#________________________________________________________________________
# Title     : get_common_column   (similar to overlay_seq_for_identical_chars )
# Usage     : %out =%{&get_common_column(\%hash1, \%hash2, '-')};
# Function  : (name1         --EHH--HHEE-- )
#             (name2         --HHH--EEEE-- ) ==> result is;
#
#             (name1_name2   -- HH--  EE-- )
#             to get the identical chars in hash strings of sequences.
#
# Example   : %out =%{&get_common_column(\%hash1, \%hash2, '-')};
#             output> with 'E' option >>> "name1     --HHH--1232-"
#   Following input will give;
#  %hash1 = ('s1', '--EHH-CHHEE----EHH--HHEE----EHH--HHEE----EHH-CHHEE--');
#  %hash2 = ('s2', '--EEH-CHHEE----EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#  %hash3 = ('s3', '-KEEH-CHHEE-XX-EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#  %hash4 = ('s4', '-TESH-CHEEE-XX-EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#
#    s1_s2_s3_s4    --E-H-CH-EE----E-H--HHEE----E-H--HHEE----E-H-CHHEE--
#
# Warning   : This gets more than 2 hashes. Not more than that!
#
# Class     : get_common_column, get_common_column_in_seq, get common column in sequence
#             for secondary structure only representation.
# Keywords  : Overlap, superpose hash, overlay identical chars, superpose_seq_hash
#             get_common_column, get_com_column
# Options   :
# Package   : Array_Util
# Reference :
# Returns   : one hash ref. of the combined key name (i.e., name1_name2). Combined by '_'
# Tips      :
# Argument  : 2 or more ref for hash of identical keys and value length.
#             One optional arg for replacing space char to the given one.
# Todo      :
# Author    : A Biomatic
# Version   : 1.5
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_common_column{
  my($i, $k,$j, $name1, $name2, @in, %out, @out_chars, $gap_chr, @str1, @str2);
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  Sub argument handling     $gap_chr here can be 'HE' etc.
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for($k=0; $k< @_ ;$k++){
	 if( ( !ref($_[$k]) )&&($_[$k]=~ /^(.)$/) ){
		$gap_chr  .= $1;    }
	 elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(.)$/) ){
		$gap_chr  .= $1;    }
	 elsif(ref($_[$k]) eq "HASH") { push(@in,  $_[$k]); }    }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"  Checking if the input hashes were right
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( (@in < 2) && ( ref($in[0]) eq 'HASH') ){
	 print "\n", __LINE__, " # get_common_column usually needs 2 hashes. Error \n";
	 print "\n", __LINE__, " # get_common_column : Anyway, I will just return the single input hash:  @in. \n";
	 %out=%{$in[0]}; # <--- This is essential to return the single input hash!!
	 goto FINAL;
  }

  %out = %{$in[0]};  ## Initializing %out
  print "\n",__LINE__, " # get_common_column hashes given are: @in \n" if $debug eq 1;

  for( $k=1; $k < @in; $k++){
	  my(@out_chars);   ## <-- Necessary to prevent concatenation.
	  my(%hash1)=%out;
	  my(%hash2)=%{$in[$k]};
	  my(@names1)= sort keys %hash1;
	  my(@names2)= sort keys %hash2;
	  $name1 = $names1[0];
	  $name2 = $names2[0];
	  @str1=split(/||\,/, $hash1{$names1[0]});
	  @str2=split(/||\,/, $hash2{$names2[0]});
	  for($i=0; $i < @str1; $i++){
	 if($str1[$i] eq $str2[$i] ){
	    push(@out_chars, $str1[$i]); }
	 elsif( defined($gap_chr) ){ push(@out_chars, $gap_chr); }
	 else{ push(@out_chars, ' '); }
	  }
	  if( $name2 < $name1){      ## To make an ordered name output eg.  seq1_seq2, than  seq2_seq1
		 %out='';
	 $out{"$name2\_$name1"}= join("", @out_chars); }
	  else{
		 %out='';
		 $out{"$name1\_$name2"}= join("", @out_chars); }
  }
  FINAL:
  \%out;
}

#________________________________________________________________________
# Title     : open_sst_files  (but reads jp file as an input!!!)
# Usage     : %out_sst_hash =%{&open_sst_files(\$jp_file_dir, \$file_name)};
# Function  : gets the name of a sst file (6taa.sst) and/or its absolute
#             dir path (aa).
#             reads the sequence names in the jp file and looks up all
#             the sst files in the same directory. Puts sst sequences
#             in a hash with keys of sequence names.
#
# Example   : jp file  ==  seq1 ABDSF--DSFSDFS   <- true sequence
#                          seq2 lkdf-jlsjlsjf
#
#                 sst files == seq1.sst, seq2.sst
#
#                 output hash == seq1 hHHHHHHHttEEEEEEEE
#                                seq2 hHHHHHHHHHEEEEEEhh
#
# Warning   : $jp_file_dir_and_name should be absolute dir and file name
#             >> This gets JP file not SST file as input !!!!
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. for a hash
# Tips      :
# Argument  : a ref. for scaler of "jp file name"
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_sst_files{
  my($sst_seqs, %jp_file, $jp, $dir,@keys,$base,$sst);
  for ($i=0; $i<=$#_; $i++){
	  if (${$_[$i]} =~ /([\w\-]+)\.sst$/){
		 $seq_name=$1;
	 ${"sst_file$i"}=${&find_seq_files($_[$i])};
	 $sst=${"sst_file$i"};
	 open(SST, "$sst");
	 while(<SST>){
	   if(/^  summary  (.*)  summary  $/){
	      $sst_seqs.=$1;	 # this is faster than    sub.
	   }
	 }
	 close (SST);
	 $out_sst_seq_hash{$seq_name }=$sst_seqs;
	  }
	  elsif (${$_[$i]} =~ /^[\w\-]+$/){
		 $dir = ${&dir_search_single($_[$i])};
		 $jp="$dir\/${$_[$i]}\.jp";
	 %jp_file=%{&open_jp_files(\$jp)};
	 @keys = (keys %jp_file);
		 for $seq_name (@keys){
	   $sst_file_name="$dir\/$seq_name\.sst";
	   open(SST, "$sst_file_name");
	   while(<SST>){
	     if(/^  summary  (.*)  summary  $/){
	       $sst_seqs.=$1;
	     }
	   }
	   close (SST);
	   $out_sst_seq_hash{$seq_name}=$sst_seqs;
	 }
	  }
	  elsif( ${$_[$i]} =~ /\/([\w\-]+\.([\w\-]+))$/){
	%jp_file=%{&read_any_seq_files($_[$i])};
	@keys = (keys %jp_file);
		for $seq_name (@keys){
	   $sst_file_name="$dir\/$seq_name\.sst";
	   open(SST, "$sst_file_name");
	   while(<SST>){
	     if(/^  summary  (.*)  summary  $/){
	        $sst_seqs.=$1;
	     }
	   }
	   close (SST);
	   $out_sst_seq_hash{$seq_name}=$sst_seqs;
	}
	  }
	  elsif( ${$_[$i]}=~ /^([\w\-]+)\.[\w\-]+$/){
	$base = $1;
		$dir=${&dir_search_single(\$base)};
	%jp_file=%{&read_any_seq_files($_[$i])};
	@keys = (keys %jp_file);
		for $seq_name (@keys){
	   $sst_file_name="$dir\/$seq_name\.sst";
	   open(SST, "$sst_file_name");
	   while(<SST>){
	     if(/^  summary  (.*)  summary  $/){
	       $sst_seqs.=$1;
	     }
	   }
	   close (SST);
	   $out_sst_seq_hash{$seq_name}=$sst_seqs;
	}
	  }
  }
  return(\%out_sst_seq_hash);
}

#________________________________________________________________________
# Title     : dir_search_single  (refer dir_search for a list of possible dirs)
# Usage     : $output_best_possible_dir = ${&dir_search_single(\$input_name)};
# Function  : With given full path or single name for a dir. it returns
#             the full path dir name. If it fails to find in pwd or given
#             specified path, it tries to search PATH, HOME etc..
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one Ref. for an array.
# Tips      :
# Argument  : One Ref. for a scalar.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub dir_search_single{    my($in_dir)=${$_[0]};
  my(@ENV_dir,$pwd,@temp,@probable_dir_list, @dirs,@possible_dirs,$final_dir_found);
  if (!(-d $in_dir)){
	  if ($in_dir =~/^[\w\.\-]+$/){ # if it is not a full path.
	  $pwd=`pwd`; chomp($pwd);
		  if(-d "$pwd\/$in_dir" ){
			 $final_dir_found = "$pwd\/$in_dir";
	  }
	  elsif( !(-d "$pwd/$in_dir" ) ){
	     @temp=split('/', $pwd); # goes up one level.
			 pop(@temp);
	     $up_pwd=join("/", @temp);
	     if (-d "$up_pwd\/$in_dir"){
	       $final_dir_found=$in_dir ;
	     }elsif( !(-d "$up_pwd/$in_dir" ) ){
				@probable_dir_list=('JPO','ALIGN','PDB','PATH','HOME','PIRDIR',
	            'PDBSST','PDBENT', 'BLASTDB', 'PIRDIR', 'SWDIR');
				for $elem (@probable_dir_list){
		   @dirs=split(/:/, $ENV{$elem});
		   for (@dirs){
		     if (/$in_dir$/){  # if $in_dir matches with a set dir.
	                $final_dir_found=$_;
	             }elsif( -d "$_\/$in_dir"){
		        $final_dir_found="$_\/$in_dir";
		     }
		   }
	        }
			 }#<---}elsif( !(-d "$up_pwd/$in_dir" ) ){
	  }
	  }elsif($in_dir =~ /\/([\w\.\-]+)$/){ # if it is a full path.
		   $in_dir = $1;
		   if(-d "$pwd\/$in_dir" ){
			  $final_dir_found = "$pwd\/$in_dir"; last;
	   }elsif( !(-d "$pwd\/$in_dir" ) ){
	      $in_dir="$up_pwd\/$in_dir";
	      $final_dir_found=$in_dir if (-d $in_dir); last
	   }else{
			  for (@probable_dir_list){
	         @dirs=split(':', $ENV{$_});
	         for (@dirs){
	           if (/$in_dir$/){
	             $final_dir_found=$_; last;
	           }
	         }
			  }
	   }
	   }
   }else{  # If the input dir is there as it is!! (no need to process)
	  $final_dir_found = $in_dir;
   }
   return(\$final_dir_found);
}
#________________________________________________________________________
# Title     : open_fas_files   (the sequence should be mininum 9 residues)
# Usage     : %anyhash = &read_fasta_files($any_sequence_file_fasta_form);
# Function  : open fasta files and put sequences in a hash
#             FASTA sequence file format is like this;
#
#             > 1st-seq
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFG
#             > 2nd.sequ
#             ABCDEFGHIJKLMOYYUIUUIUIYIKLMOPABCDEFGHIJKLMOPABCDEFG
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPPABCDEFGHIJKLMOP
#
# Example   :
# Warning   : well tested. It skips lines starting with blank, lines with '-' in them.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_fas_files{    my($input_file)=${$_[0]}; my(@names, %sequence, $flag,$name,$temp2);
   unless (-e $input_file){     print chr(7);
	 print "\n\n\t This is sub open_fas_files in $0  \n\n";
	 print "\t Fatal: The input file $input_file is not in the directory \n";  }
   open(FILE_1,"$input_file");
   while(<FILE_1>){         		# file1 needs to be xxxx.fasta for the moment, automatic later
	 if(/^\>([\w\-\.]+)\s+$/){ $name=$1; next; }
	 elsif(/^(\w\w\w\w\w\w\w\w+)$/){ $sequence{$name}.= $1; }
	 else{  next; }    }    close(FILE_1);    return(\%sequence);
}
#________________________________________________________________________
# Title     : convert_num_to_0_or_1_hash (opposite of convert_num_to_0_or_1_hash)
# Usage     : with a variable for threshold ->
#
#             %out = %{&convert_num_to_0_or_1_hash(\%input_hash, \$threshold, \%input_hash2..)};
#
# Function  : changes all the numbers into 0 or 1 according to threshold given.
#             convert_num_0_or_1_hash converts threshold and bigger nums. to
#             '0' while convert_num_0_or_1_hash_opposite converts to '1'.
# Example   : A hash =>  name1  10012924729874924792742749748374297
#                        name2  10012924729874924792710012924729874
#             A threshold => 4
#             !! if numbers are smaller than 4, they become 1 (or true).
#             Outputhash  =>  name1  11111011011111011111011011110101111
#                        name2  11111011010001011001011010010101100
#
#             ($ref1, $ref2)=&convert_num_to_0_or_1_hash(\%hash, \%hash, \$threshold);
#             above is the example when with more than 2 input hashes.
# Warning   : Threshold value is set to 0 as well as all values smaller than that.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : two references, one for hash one for scaler for threshold
#
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub convert_num_to_0_or_1_hash{
   my(@output_hash_refs, %input, $c, $i,
	  @string, $name, @names, $threshold, %output_hash);
  for($c=0; $c < @_; $c++){
	if(ref($_[$c]) eq 'SCALAR'){ $threshold =${$_[$c]};} }
  for($i=0; $i < @_; $i++){
	if(ref($_[$i]) eq 'HASH'){ %input=%{$_[$i]};
	  @names=keys (%input);
	  my($split_char)=',';
	  if ((@_ < 1)&&(ref($_[$a]) eq 'HASH')){  # if input argument is only one (= if no threshold given),
		$threshold = 1; } # <---- put 1 to $threshold as a default
	  for $name (@names){
		if($input{$name}=~/\,/){  $split_char = ','; }else{ $split_char = ''; }
		if ($input{$name} =~ /[\.\-\d]+/){ @string=split(/$split_char/, $input{$name});
	  for (@string){
			if(/\d+/){
	      if($_ > $threshold){ $_=0; } # !! becomes 0 (or false)
			  else{  $_=1;               } # !! becomes 1 (or true)
			}
		  }
		}
		$output_hash{$name}=join("", @string);
	  }
	  push(@output_hash_refs, \%output_hash);
	}
  }
  if(@output_hash_refs == 1){return($output_hash_refs[0]); }
  elsif(@output_hash_refs > 1){ return(@output_hash_refs) }
}

#________________________________________________________________________
# Title     : convert_num_0_or_1_hash_opposite (opposite of convert_num_to_0_or_1_hash)
# Usage     : with a variable for threshold ->
#
#               %out = %{&convert_num_0_or_1_hash_opposite(\%input_hash, \$threshold)};
#
# Function  : changes all the numbers into 0 or 1 according to threshold given.
#             convert_num_0_or_1_hash converts threshold and bigger nums. to
#             '0' while convert_num_0_or_1_hash_opposite converts to '1'.
# Example   : A hash =>  name1  10012924729874924792742749748374297
#                        name2  10012924729874924792710012924729874
#             A threshold => 4
#             !! if numbers are smaller than 4, they become 1 (or true).
#             Outputhash  =>  name1  11111011011111011111011011110101111
#                        name2  11111011010001011001011010010101100
#
#             ($ref1, $ref2)=&convert_num_to_0_or_1_hash(\%hash, \%hash, \$threshold);
#             above is the example when with more than 2 input hashes.
# Warning   : Threshold value is set to 0 as well as all values smaller than that.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : two references, one for hash one for scaler for threshold
#
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub convert_num_0_or_1_hash_opposite{
  my(@output_hash_refs, %input,$c, $i, $split_char,
	 @string, $name, @names, $threshold,%output_hash);
  for($c=0; $c < @_; $c++){
	if(ref($_[$c]) eq 'SCALAR'){ $threshold =${$_[$c]};} }
  for($i=0; $i <=$#_; $i++){
	if(ref($_[$i]) eq 'HASH'){
	  %input=%{$_[$i]};
	  #show_hash(\%input);
	  @names=keys (%input);
	  $split_char=',';
	  if ((@_ == 1)&&(ref($_[$a]) eq 'HASH')){  # if input argument is only one (= if no threshold given),
		$threshold = 1; } # <---- put 1 to $threshold as a default
	  for $name (@names){
		if($input{$name}=~/\,/){  $split_char = ',';
		}else{ $split_char = ','; }
		if ($input{$name} =~ /[\.\-\d]+/){
		  @string=split(/$split_char/, $input{$name});
	  for (@string){
			if(/\d+/){
	      if($_ >= $threshold){  # !! becomes 1 (or false)
				 $_ = 1;
			  }else{  $_=0;   } # !! becomes 0 (or true)
			}
		  }
		}
		$output_hash{$name}=join(",", @string);
	  }
	  push(@output_hash_refs, \%output_hash);
	}
  }
  #show_hash(\%output_hash);
  if(@output_hash_refs == 1){return($output_hash_refs[0]); }
  elsif(@output_hash_refs > 1){ return(@output_hash_refs) }
}

#________________________________________________________________________
# Title     : superpose_seq_hash   ( first to second hash)
# Usage     : %out =%{&superpose_seq_hash(\%hash1, \%hash2)};
# Function  : (name1 000000112324)+(name1  ABC..AD..EFD ) => (name1 000..01..324)
#             uses the second hash a template for the first sequences. gap_char is
#             '-' or '.'
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one hash ref.
# Tips      :
# Argument  : 2 ref for hash of identical keys and value length.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub superpose_seq_hash{ my($gap_chr)='-'; my($i,$j,$name,%out,@str1,@str2);
  if((ref($_[0]) eq 'HASH')&&(ref($_[1]) eq 'HASH')){
	my(%hash1)=%{$_[0]}; my(%hash2)=%{$_[1]}; }
  else{ print "\n superpose_seq_hash needs hash ref\n"; print chr(007); exit; }
  my(@names1)=keys %hash1;my(@names2)=keys %hash2;
  (@names1 > @names2)? $bigger=@names1 : $bigger=@names2;
  for ($j=0; $j < $bigger; $j++){
	if($hash2{$names2[$j]}=~/(\W)/){ $gap_chr = $1; }
	  @str1=split(//, $hash1{$names1[$j]}); @str2=split(//, $hash2{$names2[$j]});
	  for($i=0; $i < @str2; $i++){
		if(($str2[$i] =~ /\W/)||($str2[$i] =~ //)){ $str1[$i]="$gap_chr";}}
	  $out{$names1[$j]}=join(",", @str1);
  }
  \%out;
}

#________________________________________________________________________
# Title     : open_fasta_files   (the sequence should be mininum 9 residues)
# Usage     : %anyhash = &read_fasta_files($any_sequence_file_fasta_form);
# Function  : open fasta files and put sequences in a hash
#             FASTA sequence file format is like this;
#
#             > 1st-seq
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFG
#             > 2nd.sequ
#             ABCDEFGHIJKLMOYYUIUUIUIYIKLMOPABCDEFGHIJKLMOPABCDEFG
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPPABCDEFGHIJKLMOP
#
# Example   :
# Warning   : well tested. It skips lines starting with blank, lines with '-' in them.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_fasta_files{    my($input_file)=${$_[0]}; my(@names, %sequence, $flag,$name,$temp2);
   unless (-e $input_file){     print chr(7);
	 print "\n\n\t This is sub open_fas_files in $0  \n\n";
	 print "\t Fatal: The input file $input_file is not in the directory \n";  }
   open(FILE_1,"$input_file");
   while(<FILE_1>){         		# file1 needs to be xxxx.fasta for the moment, automatic later
	 if(/^\>([\w\-\.]+)\s+$/){ $name=$1; next; }
	 elsif(/^(\w\w\w\w\w\w\w\w+)$/){ $sequence{$name}.= $1; }
	 else{  next; }    }    close(FILE_1);    return(\%sequence);
}
#________________________________________________________________________
# Title     : open_pdb_files  (read the sequences only)
# Usage     : %out = %{&open_pdb_files(\$VAR)};
# Function  : Convert a PDB structure file to FASTA format sequences.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : One ref. for a hash of sequences(DNA, RNA, PROTEIN (IN diff chains)
#             If the two chains are identical, it rids of one of them and returns
#             a name with out chain note-->  2tma, not 2tmaA and 2tmaB
# Tips      :
# Argument  : one ref. for an inputfile (absolute
#             >>> PDB example >>>
#             SEQRES   1 A  284  MET ASP ALA ILE LYS LYS LYS MET GLN MET LEU LYS LEU  2TMA  51
#             SEQRES   2 A  284  ASP LYS GLU ASN ALA LEU ASP ARG ALA GLU GLN ALA GLU  2TMA  52
#
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_pdb_files{    my($input)=${$_[0]};  $input=${&find_seq_files(\$input)};
  my($outseq, $structure, %outhash, @fields, %AA, $chain, $residues);
	$AA{"ALA"} = "A";  $AA{"MET"} = "M";  $AA{"ASP"} = "D";  $AA{"PRO"} = "P";
	$AA{"CYS"} = "C";  $AA{"ASN"} = "N";  $AA{"GLU"} = "E";  $AA{"GLN"} = "Q";
	$AA{"PHE"} = "F";  $AA{"ARG"} = "R";  $AA{"GLY"} = "G";  $AA{"SER"} = "S";
	$AA{"HIS"} = "H";  $AA{"THR"} = "T";  $AA{"ILE"} = "I";  $AA{"VAL"} = "V";
	$AA{"LYS"} = "K";  $AA{"TRP"} = "W";  $AA{"LEU"} = "L";  $AA{"TYR"} = "Y";
  open(INPUT_PDB_FILE, "$input");
  while (<INPUT_PDB_FILE>){
	if(/^HELIX/ || /^ATOM/ || /^FTNOTE/ || /^HET/){ last; };
	if(/^SEQRES +\d+ +(.) +(\d+) +(.+)\s+(\w+)\s+\d+.+$/){
	  $chain =$1; $seq_size = $2; $residues = $3; $pdb_name=$4;
	  @residues=split(' ', $residues);
	  if($residues[0]=~/[A-U]/){ # <-- Check if it is DNA/RNA seq.
		for(@residues){  $outhash{"\L$pdb_name\U$chain\E"}.=$_;  }
	  }else{
	for (@residues){  $outhash{"\L$pdb_name\U$chain\E"}.=$AA{$_}; } } }
  }
  close( INPUT_PDB_FILE);
  @keys=keys %outhash;
  if ($chain=~/[A-Z]/){
	for ($i=0; $i < @keys; $i++){
	  for ($j=$i+1; $j < @keys; $j++){
		if ($outhash{$keys[$i]} eq $outhash{$keys[$j]}){
		  delete($outhash{$keys[$j]});  $temp=$keys[$i];
	  chop($temp);  $outhash{$temp}=$outhash{$keys[$i]};
	  delete($outhash{$keys[$i]});	} } } }  return( \%outhash  );
}
#________________________________________________________________________
# Title     : open_msf_files
# Usage     : (*out, *out2) = @{&open_msf_files(\$inputfile1, \$inputfile2)};
#             : (@out)        = @{&open_msf_files(\$inputfile1, \$inputfile2)};
#             ---------- Example of MSF ---
#             PileUp
#
#             MSF:   85  Type: P    Check:  5063   ..
#
# Function  : open msf files and put sequences in a hash(s)
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : (*out, *out2)  or (@out_array_of_refs)
# Tips      :
# Argument  : (\$inputfile1, \$inputfile2, .... )};
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_msf_files{  my(@in)=@_; my(@names, $n, $s, %hash,@out_hash_ref_list);
  for($i=0; $i<=$#in; $i++){
	if(ref($in[$i])){ unless (-e ${$_[$i]}){ next; }
	   open(FILE_1,"${$_[$i]}");  undef(%hash);
	   while(<FILE_1>){      # file1 needs to be xxxx.msf for the moment, automatic later
	  if((/^([\S]+)\t* +$/)||(/^\#/)||(/^\-+/)){ next; }
	  if(/^([\S]+)\t* +([\.\w ]+)[\n]$/){ $n=$1;$s=$2; $s=~s/ //g; $hash{$n}.= $s; }}
	   close(FILE_1);  push(@out_hash_ref_list, \%hash); } }
  if($#out_hash_ref_list  == 0 ){ return(\%hash); }
  elsif($#out_hash_ref_list > 0){ return(@out_hash_ref_list); } # <-- contains (\%out_seq0, \%out_seq1, \%out_seq2, .... )
}
#________________________________________________________________________
# Title     : hash_output_chk
# Usage     : &hash_output_chk(\%outing_hash);
# Function  : checks hash output of any subroutine.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub hash_output_chk{
  for($i=0; $i<= $#_; $i++){
	my(%tem)=%{$_[$i]};  my(@keys)=keys %tem;
	for ($j =0; $j<@keys; $j++){
	  unless(($keys[$j]=~/[\s\S]+/)&&($tem{$keys[$j]}=~/[\s\S]+/)){
		print "\n Err. at Hash_output_chk at $0 \n", chr(7); exit;
	  }
	}
  }
}

#________________________________________________________________________
# Title     : open_slx_files
#             selex file (foo.slx) looks like this:
#
#             #=SQ GLB_TUBTU  5.9393 - - 0..0::0 -
#             #=SQ GGZLB      20.9706 - - 0..0::0 -
#
#             HAHU        .........VLSPADKTNVKAAWGKVGA......HAGEYGAEALERMFLS
#             HBA3_PANTR  .........VLSPADKTNVKAAWGKVGA......HAGZYGAEALERMFLS
#
# Usage     : %anyarray = &open_slx_files(\$any_sequence_file_msf_form);
# Function  : open slx files and put sequences in a hash
# Example   :
# Warning   : The MSF FORMAT SHOULD BE AT LEAST 30 residue long
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. of a hash
# Tips      :
# Argument  : takes one ref. for a file.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_slx_files{    	                    my(@names, $n, $s, %hash);
   if ((-z ${$_[0]})|| (-B _) || (-x _)){   print chr(7);
	  print "\n\t I am $0: Input file $file1 isn't in the dir \n"; exit;
   }
   open(FILE_1,"${$_[0]}");  	# reading in (MSF)
   while(<FILE_1>){         	# file1 needs to be xxxx.msf for the moment, automatic later
	 if((/^([\S]+)[\t]* +$/)||(/\#/)){ next; }
	 if(/^([\w\_\.\=\+\#\@]+)\t* \t*([\.\w]+)[\n]$/){
		$n=$1; $s=$2; $hash{$n}.= $s;  }
   }close(FILE_1);  return( \%hash );
}
#________________________________________________________________________
# Title     : open_jp_files
# Usage     : %out_hash=%{&open_jp_files(\$file_name)};
# Function  : reads jp files and stores results in a hash.
# Example   :
# Warning   : All the spaces  '-' !!!
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a reference of a hash for names and  their sequences.
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_jp_files{      my(%hash_out, $s, $n, $s1);
   open(FILE_JP, "${$_[0]}");  # reading in (JP) file
   while(<FILE_JP>){    if(/^CLUSTAL/){ next; } # <-- a safety measure
	 if((/^([\S]+)[\t]* +$/)||(/\#/)){ next; }
	 if(/^([\w\.\-\=\+]+) +\t*(.+)[\n]$/){ $n=$1; $s=$2; $hash_out{$n}.= $s; }
   }close(FILE_JP);  \%hash_out;
}

#________________________________________________________________________
# Title     : open_aln_files
# Usage     : %out_hash=%{&open_jp_files(\$file_name)};
# Function  : reads jp files and stores results in a hash.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a reference of a hash for names and  their sequences.
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_aln_files{       my(%hash_out, $s, $n, $s1);
  open(FILE_JP, "${$_[0]}");  # reading in (JP) file
  while(<FILE_JP>){       if(/^CLUSTAL/){ next; }
	if((/^([\S]+)\t* +$/)||(/\#/)){ next; }
	if(/^([\w\.\-\=\+]+) +\t*(.+)[\n]$/){ $n=$1; $s=$2; $hash_out{$n}.= $s; }
  }close(FILE_JP);   \%hash_out;
}
#________________________________________________________________________
# Title     : get_base_name
# Usage     : $dir = ${&pwd_dir($any_absolute_path_dir)};
# Function  : returns present working dir base
# Example   :
# Warning   : well tested.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_base_name{
  my(@pwd)=split(/\//,`pwd`);
  my($dir)=$pwd[$#pwd];
  chomp($dir);	# !!! chomp() is essential, as `pwd` returns new line.
  \$dir;          # this should be present unless chomp returns 1 !!!
}
#________________________________________________________________________
# Title     : open_sst_files_with_gap  (but reads jp file as an input, too!!!)
# Usage     : %out_sst_hash =%{&open_sst_files_with_gap(\$jp_file_dir_and_name)};
# Function  : gets the name of a file(jp file) with its absolute dir path
#             reads the sequence names in the jp file and looks up all
#             the sst files in the same directory. Puts sst sequences
#             in a hash with keys of sequence names.
#
# Example   : jp file  ==  seq1 ABDSF--DSFSDFS   <- true sequence
#                              seq2 T--kdf-GAGGGASF     (aligned)
#
#                 sst files ==> 'seq1.sst', 'seq2.sst' (in the same dir)
#
#             original sst format:  seq1 hHHHHHttEEEE  <-- No gaps!
#                                  seq2 hHHHHHHEEhh
#             After this sub ==>
#             (final out hash =   (  seq1 hHHHH--HttEEEE  <-- inserted
#                                  seq2 h--HHH-HHHEEEhh  )     gaps
#
# Warning   : $jp_file_dir_and_name should be absolute dir and file name
#             >> This gets JP file not SST file as input !!!!
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. for a hash
# Tips      :
# Argument  : a ref. for scaler of "jp file name"
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_sst_files_with_gap{   # This automatically determines MSF or JP format
  my(%seq_file, $sst_file, %secondary_struc, @keys, $directory);
  %seq_file = %{read_any_seq_files($_[0])};

  ######  Simple reading in of SST files ------------
  ######  Simple reading in of SST files ------------
  for $seq_name (keys %seq_file){
	 $sst_file ="$seq_name\.sst";
	 print $sst_file;
	 %secondary_struc =( %secondary_struc, %{&read_any_seq_files(\$sst_file)} );
  }
  print %secondary_struc;
  ### Now we have  1. %jp_file  and  2. %out_sst_seq_hash  -------
  if (!(defined(%secondary_struc))){
	 return(\%seq_file);
  }else{
	%gap_corrected_out=%{&put_gaps_in_hash(\%seq_file, \%secondary_struc)};
	return( \%gap_corrected_out );
  }
}
#________________________________________________________________________
# Title     : put_gaps_in_hash  (The order of input hashes does matter)
# Usage     : %out=%{&put_gaps_in_hash(\%hash_with_gap, \%hash_sans_gap)};
#
#             %hash1=('1ctx',  '111111111111111',      <-- hash input without gaps
#             '2ctx',  '2222222222222222',
#             '3ctx',  '3333333333');
#
#             %hash2=('1ctx',  'AAA--AAAAAAAAAAAA-',   <-- hash input with template gaps
#             '2ctx',  'BBBBBBBBBBBB-BBBB',
#             '3ctx',  'CCCCCC----CCCC');
#
#             >> resulting out hash;
#
#             %hash3=('1ctx',     '111--111111111111-',
#             '2ctx',     '222222222222-2222',
#             '3ctx',     '333333----3333 );
#
# Function  :
# Example   :
# Warning   : The keys for hashes should be the same and the two sequences
#             should be identical.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one hash reference.
# Tips      :
# Argument  : 2 hash references.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub put_gaps_in_hash{
  my($temp0)= values %{$_[0]};  # finds the hash with gaps
  my($temp1)= values %{$_[1]};  # and assigns to right input hash.
								# above puts the first values to $temp0 & 1
  if (($temp0=~/\-/)||($temp0=~/\./)){  # compares the leng of the first
	%hash_gap = %{$_[0]};               # values of hashes and assigns
	 %hash_sans_gap=%{$_[1]};            # accordingly.
  }elsif(($temp1=~/\-/)||($temp1=~/\./)){
	%hash_gap     =%{$_[1]};
	%hash_sans_gap=%{$_[0]};
  }else{
	%hash_gap     =%{$_[0]};  # If it can not determine input type, it assumes
	%hash_sans_gap=%{$_[1]};  # that the first one was for gap, the 2nd for secondary.
  }                           # structure or whatever.
  my(@keys)=sort keys (%hash_gap);
  my($gap_char) = '-';  #  default gap_char is  '-'
  my(@string1, @string2, @gap_pos, %out_hash, $gapped_string, $res);

  if ($hash_gap{$keys[0]}=~/\-/){
	$gap_char = '-';
  }elsif($hash_gap{$keys[0]}=~/\./){
	$gap_char = '.';
  }
  ########## Actual exchange part ############
  for (@keys){
	@string1 = split('', $hash_gap{$_});
	@string2 = split('', $hash_sans_gap{$_});
	 for ($t=0; $t <=$#string1; $t++){
	   $res=$string1[$t];
		if(($res =~ /\-/)||($res =~ /\./)||($res =~ /\s/)){
		  splice(@string2, $t, 0, $gap_char);
		}
	 }
	$gapped_string = join("", @string2);
	$out_hash{$_}= $gapped_string;
  }
  \%out_hash;
}
#________________________________________________________________________
# Title     : hash_common_by_keys
# Usage     : %hash1_value = %{&hash_common_by_keys(\%hash1, \%hash2,...)};
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : the VALUES OF THE FIRST HASH which occur in later hashes
#             are returned
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub hash_common_by_keys{  my(%common)=();
  for($i=0; $i< @_; $i++){  my(%common2)=();
	 if( !(defined(%common) )){ %common=%{$_[$i]}; next;}
	 elsif(defined(%common)){ %h1=%{$_[$i]};
	   for(keys %common){ $common2{$_}=$common{$_} if (defined $h1{$_});}
	 %common=%common2;}
	 undef(%common2);  }
  \%common;
}

#________________________________________________________________________
# Title     : get_base_names
# Usage     : $base =${&get_base_names(\$file_name)};
#             :   or @bases = &get_base_names(\@files);  # <-- uses `pwd` for abs directory
# Function  : produces the file base name(eg, "evalign"  out of "evalign.pl" ).
# Example   : $base => 'test'  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Class     :
# Keywords  : get_base_name{, base_name, file_base_name ,  get_file_base_name
# Options   :
# Package   : Jong_Util
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_base_names{
	my($x, @out_file, $file, @file, $base, @base);
	@file=@{$_[0]} if (ref($_[0]) eq 'ARRAY');
	@file=@_ if !(ref($_[0]) eq 'ARRAY');
	for($x=0; $x < @file; $x ++){
		if( ref($file[$x]) ){
			$file = ${$file[$x]};
			@out_file=split(/\./, $file);
			$base = $out_file[0];
		}else{
			$file = $file[$x];
			@out_file=split(/\./, $file);
			$base = $out_file[0];
		}
		push(@base, $base);
	}
	if(@base == 1 ){ \$base[0] }else{ @base }
}

#________________________________________________________________________
# Title     : get_posi_sans_gaps (get positions without after removing gaps)
# Usage     : @seq_position1 = &get_posi_sans_gaps($string1);
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : the positions of residues after removing gaps(but keeps pos).
#               used for analysis of shifted positions of seq. comparison.
# Tips      :
# Argument  : one scalar variable input of sequence string.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_sans_gaps{
  my($string_input) = ${$_[0]};
  my($char, @positions, $i);
  for($i=0; $i < length($string_input); $i++){
	$char=substr($string_input,$i,1);
	if(($char eq '-')||($char eq '.')){  next; }else{ push(@positions, $i); } }
  \@positions;
}


#________________________________________________________________________
# Title     : get_posi_diff    # used in 'get_posi_shift_hash'
# Usage     : @position_diffs =&get_posi_diff(\@seq_position1,\@seq_position2);
# Function  :
# Example   : @compacted_posi_dif =(1 ,2, 1, 1, '.' ,2,  1,  1, '.');
#             @compacted_posi_dif2=(4 ,2, 1, 1, ,2,  1, '.' ,3,  1);
#             output ==> ( 3 0 0 0 . 1 . 2 .)   (it ignores positions which have non digits.
#             output ==> (-3 0 0 0 . 1 .-2 .) when abs is not used.
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. for an @array of differences of input arrays. array context.
# Tips      :
# Argument  : Takes two ref. for arrays which have positions of residues.
# Todo      :
# Author    : A Biomatic
# Version   : 1.4
# Used in   : evalign.pl, get_position_shift_rate
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_diff{
   my(@positions1)=@{$_[0]};
   my(@positions2)=@{$_[1]};
   my(@num_diffs_between_str_and_ali, $diff, $z, $gap_char);
   if($debug eq 1){
	 print __LINE__, " # get_posi_diff : \n";
   }
   $gap_char = '.';
   for ($z=0; $z < @positions2; $z++){
	 if (($positions1[$z] =~ /\d+/) && ($positions2[$z] =~ /\d+/)){
		$diff=($positions1[$z] - $positions2[$z]);
	push(@num_diffs_between_str_and_ali, $diff );
	 }else{
		push(@num_diffs_between_str_and_ali, $gap_char);
	 }
   }
   \@num_diffs_between_str_and_ali;
}
#________________________________________________________________________
# Title     : remov_com_column
# Usage     : %new_string = %{&remov_com_column(\%hashinput)};
#             @out=@{&remov_com_column(\@array3)};
# Function  : removes common gap column in seq.
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. of  hash(es) and array(s).
#
#             name1   ABCDE
#             name2   ABCDEE
#             name3
#
#             (ABC
#             from above the two column of dot will be removed
#             To remove absurd gaps in multiple sequence alignment
# Tips      :
# Argument  : accepts reference for hash(es) and array(s).
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remov_com_column{ my(@hash_ref_out, $d, $gap_char);
  for($d=0; $d < @_; $d++){
	 if(ref($_[$d]) eq 'HASH'){ my($len,@string,@new_string,@string2);
		my(%input)=%{$_[$d]};  my(@common);
		for (keys %input){
	   @string = split('', $input{$_});
		   if(!(defined(@common))){ @common = (@string);  }
		   else{ for ($k=0; $k < @string; $k++){
			 if (($string[$k] =~ /\W/ )&&($common[$k] =~ /(\W)/)){ $common[$k]= $1;}
			 elsif(($string[$k] =~ /(\W)/)&&(!(defined($common[$k])))){ $common[$k]=$1;}
			 else{ $common[$k]='X';} } } }
		for (keys %input){ @string2 = split(//, $input{$_});
		   for ($i=0; $i < @string2; $i++){
			 if ($common[$i] ne $string2[$i]){ push(@new_string, $string2[$i]); } }
		   $new_string{$_}= join('', @new_string); @new_string = ();      }
		push(@hash_ref_out, \%new_string);
	 }
	 elsif(ref($_[$d]) eq 'ARRAY'){  my( $k, $y, $x,@string_array, @string);
		my(@input)=@{$_[$d]};  @common=();
	for($y=0; $y< @input; $y++){
	   @string = split('', $input[$y]);
		   if(!(defined(@common))){  @common = @string;  }
		   else{
	     for ($k=0; $k < @string; $k++){
				if (($string[$k]  =~ /(\W)/)&&($common[$k]  =~ /(\W)/)){ $common[$k]=$1;}
				elsif(($string[$k] =~ /(\W)/)&&(!(defined($common[$k])))){ $common[$k]=$1;}
				else{ $common[$k]='X';}
	     }
	   }
	}
	for($x=0; $x < @input; $x++){  my($new_string, @string2);
	   @string2 = split(//, $input[$x]);
		   for ($i=0; $i < @string2; $i++){
			  if ($common[$i] ne $string2[$i]){ $new_string.= "$string2[$i]"; }
		   }
	   push(@string_array, $new_string);
	}
		push(@hash_ref_out, \@string_array);
	 }
  }
  if(@hash_ref_out ==1) { return( $hash_ref_out[0] ); }
  elsif(@hash_ref_out>1){ return( @hash_ref_out ) }
}
#________________________________________________________________________
# Title     : put_position_back_to_str_seq ( put_posi_back_to_str_seq )
# Usage     : @result =@{&put_position_back_to_str_seq(\@string_from_struct, \@compacted_posi_dif)};
# Function  :
# Example   : @string_from_struct=('X', 'T', 'A' ,'B' , '.' ,'F',  'G', '.' , 'O' ,'P', '.');
#             @compacted_posi_dif=(1 ,2, 1, 1, ,2, 1, 1, 1);
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. for an array
# Tips      :
# Argument  : takes two refs for arrays (one for char the other for digits
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub put_position_back_to_str_seq{
  my(@string_from_struct)=@{$_[0]};
  my(@compacted_posi_dif)=@{$_[1]};
  my($j)=0; my($char)=0; my($i);
  for ($i=0; $i < @string_from_struct; $i++){
	$char = $string_from_struct[$i];
	if ($char =~ /\w/){
	   $string_from_struct[$i] = $compacted_posi_dif[$i-$j];
	}else{ $j++; }
  }
  \@string_from_struct;
}


#________________________________________________________________________
# Title     : find_seq_files (similar to find.pl) used in 'read_any_seq_file.pl'
# Usage     : $found_file = ${&find_seq_files(\$input_file_name)};
# Function  : seeks given test file in pwd, specified dir, default path etc.
#             If not found yet, it looks at all the subdirectories of path and pwd.
#             PATH environment dirs, then returns full path file name.
# Example   : $found_file=${&find_seq_files(\$input_file_name)};
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : return( \$final );
# Tips      :
# Argument  : (\$input_file_name) while $input_file_name can be  'xxx.xxx', or '/xxx/xxx/xxx/xxy.yyy'
#             or just directory name like 'aat' for  /nfs/ind4/ccpe1/people/A Biomatic /jpo/align/aat
#             then, it tries to find a file with stored seq file extensions like msf, jp, pir etc
#             to make aat.msf, aat.jp, aat.pir ... and searches for these files.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub find_seq_files{
  my($final, $no_ext_file, $result); my($in_file)=${$_[0]}; my($pwd)=`pwd`; chomp($pwd);
  my( $base, @ENV_dir, $ext, @probable_dir_list, $directory);
  my(@extension_db)=('sst','msf','fasta','jp','fas','aln','brk','pdb', 'rms', 'ent','slx','fa');
  my(@probable_dir_list)=('JPO','ALIGN','PATH','HOME','PIRDIR','PWD','PDBSST','PDBENT','BLASTDB','PIRDIR','SWDIR','PDB');
   if(($in_file=~/\//)&&(-e $in_file)){ $final=$in_file; }
   elsif((-e $in_file)&&(-s $in_file)&&($in_file !~/\//)){ $in_file="$pwd\/$in_file"; $final=$in_file;}
   ######## if it was like  '/nfs/ind4/ccpe1/people/A Biomatic /perl.msf'
   elsif($in_file =~ /\/([\w\-\.]+)$/){ $in_file = $1;
		if(-e $in_file){ $final = "$pwd\/$in_file"; }
	#### if it has xxxxxx.xxxx  file form. #######
	elsif($in_file =~ /(([\w\-]+)\.([\w\-]+))$/){ $file=$1; $base=$2; $ext=$3;
	    for (@extension_db){ if($_ eq $ext){ shift(@extension_db);}}
	    unshift(@extension_db, $ext);
	    for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
	       push( @ENV_dir, split(/:/, $ENV{$_}));}
	       for $dir (@ENV_dir){ $in_file="$dir\/$file";
		  if ((-e $in_file) && (-s $in_file)){  $final=$in_file; last;}
		  else{
		      for $ext (@extension_db){ $in_file="$dir\/$base\.$ext";
			  if ((-e $in_file) && (-s $in_file)){
			     if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
	       unless(defined ($final)){
		  for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
		      if(-e $in_file){ $final=$in_file; last; }}}}

	 ### if it has  xxxxxx   file form, ie. not extension #######
	 elsif($in_file =~/\/([\w_\-]+)$/){  $base = $1;
	   for (@extension_db){
	     if($_ eq $ext){ shift(@extension_db);  }
	     unshift(@extension_db, $ext);
	     for (@probable_dir_list){
	       if ($ENV{$_} =~ /\/$/){  chop($ENV{$_}); }
	       push( @ENV_dir, split(/:/, $ENV{$_}) );
	       for $dir (@ENV_dir){ $no_ext_file="$dir\/$base";
		   if((-e $no_ext_file) && (-s $no_ext_file)){ $final=$no_ext_file; last;}
		   else{
		     for $ext (@extension_db){ $in_file ="$dir\/$base\.$ext";
		        if ((-e $in_file) && (-s $in_file)){ $final = $in_file; last;}}}}}}}}

	#### when the input was like this  'perl.msf'  in any directory.
	elsif($in_file =~ /^(([\w\-]+)\.([\w\-]+))$/){ $file=$1; $base=$2; $ext=$3;
	for (@extension_db){ if($_ eq $ext){ shift(@extension_db);}}
	unshift(@extension_db, $ext);
	for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
	   push( @ENV_dir, split(/:/, $ENV{$_}));}
	   for $dir (@ENV_dir){ $in_file="$dir\/$file";
	      if ((-e $in_file) && (-s $in_file)){ $final=$in_file; last;}
	      else{
		  for $ext (@extension_db){ $in_file="$dir\/$base\.$ext";
		      if ((-e $in_file) && (-s $in_file)){
		         if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
	   unless(defined ($final)){
	      for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
		  if(-e $in_file){ $final=$in_file; last; }}}}
	#### when the input was like this  'hemocyan'  in any directory.
	elsif($in_file =~ /^([\w\-]+)$/){ $file=$1;
	for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
	   push( @ENV_dir, split(/:/, $ENV{$_}));}
	   for $dir (@ENV_dir){ $in_file="$dir\/$file";
	      if ((-e $in_file) && (-T $in_file)){  $final=$in_file; last;}
	      else{
		  for $ext (@extension_db){ $in_file="$dir\/$file\.$ext";
		      if ((-e $in_file) && (-s $in_file)){
		         if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
	   unless(defined ($final)){
	      for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
		  if(-e $in_file){ $final=$in_file; last; }}}}
   END_POINT:
   return( \$final );
}
#________________________________________________________________________
# Title     : read_any_seq_files  (can handle multiple input)
# Usage     : %out_seq=%{&read_any_seq_file(\$input_file_name)};
# Function  : Tries to find given input regardless it is full pathname, with or
#             without extension. If not in pwd, it searches the dirs exhaustively.
# Example   : (*out1,  *out2) =&read_any_seq_files(\$input1, \$input2);
#             : (@out_ref_array)=@{&read_any_seq_files(\$input1, \$input2)};
#             : (%one_hash_out) =%{&read_any_seq_files(\$input1)};
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : 1 ref. for a HASH of sequence ONLY if there was one hash input
#             1 array (not REF.) of references for multiple hashes.
# Tips      :
# Argument  : one of more ref. for scalar.
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_any_seq_files{
  my(@out_hash_ref_list, $o, $ext );
  my(@in)=@_;
  for($o=0; $o< @in; $o++){
	my($found, %out, @file_ext_accepted, $found_file, $sub);
	if(ref($_[$o])){
	   @file_ext_accepted=('msf', 'fasta','jp','aln','ali','pir',
						  'slx', 'dna','fas','pdb','rms','brk', 'dssp');
	   if( !(-e ${$in[$o]} )||(-B ${$in[$o]})||(-z ${$in[$o]} ) ){
		  print "\nSUB read_any_seq_files: ${$in[$o]} no seq file exists(or not passed at all) for $0 \n\n",
		  chr(7); exit;}
	   $found_file=${&find_seq_files($in[$o])};
	   for $ext(@file_ext_accepted){ my($sub)="open\_$ext\_files";
		 if($found_file =~/\.$ext$/){%out=%{&{"$sub"}(\$found_file)} if (defined &{"$sub"}); $found =1;}}
	   if($found==0){ my($sub)="open\_$ext\_files"; #<--- this is the last resort !!
		 for $ext(@file_ext_accepted){%out=(%out, %{&{"$sub"}(\$found_file)}) if (defined &{"$sub"});}}}
	   elsif(!(ref($_[$o]))){ print "\nread_any_seq_files in $0 files accepts only REFERENCES\n\n"; exit;}
	 push(@out_hash_ref_list, \%out);
  }
  if(@out_hash_ref_list == 1){  ### If only single hash output is,
	 return($out_hash_ref_list[0]);
  }elsif( @out_hash_ref_list > 1){
	 return(@out_hash_ref_list);   # <-- contains (\%out_seq0, \%out_seq1, \%out_seq2, .... )
  }
}

#________________________________________________________________________
# Title     : dir_name  (same as  pwd_dir)
# Usage     : $dir = &pwd_dir($any_absolute_path_dir);
# Function  : returns present working dir name
# Example   :
# Warning   : well tested.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub dir_name{ my(@pwd)=split(/\//,`pwd`); my($dir)=$pwd[$#pwd]; chomp($dir); \$dir;
}

#________________________________________________________________________
# Title     : get_base_name
# Usage     : $base =${&get_base_name(\$absolute_dir)};
#             or $base =${&get_base_name};  # <-- uses `pwd` for abs directory
# Function  : produces the file base name.
# Example   : $base => test  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_base_name{  my(@pwd);   if ($#_ == 0 ){ @pwd=split(/\//, ${$_[0]});}
   else{ @pwd=split(/\//,`pwd`); } my($dir)=$pwd[$#pwd]; chomp($dir);  \$dir;
}


#________________________________________________________________________
# Title     : dir_path  (same as  pwd_path )
# Usage     : $any_path = &dir_path($any_directory); or &dir_path('.') for pwd.
# Function  : returns directory path (= pwd ), eg.  /nfs/ind4/ccpe1/people/A Biomatic
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub dir_path{   my($pwd)=`pwd`;  chomp($pwd);	return( \$pwd ); }



#________________________________________________________________________
# Title     : get_each_posi_diff_hash  (used in get_posi_rates_hash_out)
# Usage     : %position_diffs =%{&get_each_posi_diff_hash(\@seq_position1, \@seq_position2)};
# Function  : gets a ref. of a hash and calculates the position diffs.
# Example   :
# Warning   : split and join char is ',';
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. for an array of differences of input arrays. array context.
#             ---Example input (a hash with sequences); The values are differences after
#                                comparion with structural and sequential alignments.
#             %diffs =('seq1', '117742433441...000',   <-- input (can be speparated by '' or ','.
#               'seq2', '12222...99999.8888',
#             'seq3', '66222...44444.8822',
#             'seq4', '12262...00666.772.');
#             example output;
#             seq3_seq4       '0,1,0,0,0,.,.,.,,.,0,,0,0,,0,0,,.,0,,0,0,.'
#             seq1_seq2       '0,1,0,1,1,.,.,.,,.,2,,2,2,,2,2,,.,.,,2,2,1'
#             seq1_seq3       '0,1,0,1,1,.,.,.,,.,1,,1,1,,0,.,,.,.,,1,1,1'
#             seq1_seq4       '0,1,0,,1,1,.,.,.,,.,1,,1,1,0,.,.,,.,1,,2,2'
#             seq2_seq3       '0,1,0,,0,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,0'
#             seq2_seq4       '0,0,0,,1,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,.'
# Tips      :
# Argument  : Takes a ref. for hash which have positions of residues of sequences.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_each_posi_diff_hash{
   my(%diffs)= %{$_[0]};
   my(@names)= keys (%diffs);
   my(%seqs_comp_in_pair, $i, $j, $k, @temp, @temp2,$split_char);
   for ($i=0; $i < @names; $i++){
	  if($diffs{$names[$i]}=~/\,/){ $split_char =',';}else{ $split_char = ''; }
	  (@{"string$i"}) = split(/$split_char/, $diffs{$names[$i]});   }
   for ($i=0; $i < @names; $i++){
	  for ($j=$i+1; $j < @names; $j ++){
	 for ($k=0; $k < @string0; $k++){
			if ((${"string$i"}[$k] =~ /[-\d+]/) && (${"string$j"}[$k] =~ /[-\d+]/)){
	       my($diff) = abs(${"string$i"}[$k] - ${"string$j"}[$k]);
	       push(@temp2, $diff); }
	    else{ push(@temp2, '.'); } }

		 ####################################################################################
		 ###### Following if {} is for sorting output names to make  2aaa_6taa than 6taa_2aaa
		 ####################################################################################
		 if($names[$i] <= $names[$j]){
			$seqs_comp_in_pair{"$names[$i]\_$names[$j]"}=join(",", @temp2); }
		 else{ $seqs_comp_in_pair{"$names[$j]\_$names[$i]"}=join(",", @temp2); }

		 @temp2=();
	  }
	}
   \%seqs_comp_in_pair;  # permutated
}

#________________________________________________________________________
# Title     : sum_hash
# Usage     : $out = &sum_array(%anyhash);
# Function  : sum of all the  numbers in valuse of a hash
# Example   : %hashinput= ( name1, '12..3e',
#                            name2, '...234');
#             $result = 1+2+3+2+3+4 = 15 (from above example)
# Warning   : It only gets digits in the input strings and sums them up.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. of a scaler.
# Tips      :
# Argument  : ref. of an array of numbers.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub sum_hash{		# usage:  $output = &sum_array(@any_array);
	my($elements) = join(",", values (%{$_[0]}));
	my(@elements) = split(',',$elements);
	my($sum);
	foreach $item(@elements){
	  if ($item =~ /[\-\d+]/){
		 $sum += $item;
	  }
	}
	\$sum;
}
#________________________________________________________________________
# Title     : sum_array
# Usage     : $out =  ${&sum_array(\@anyarray)};
# Function  : sum of all the  elements of an array .
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. of a scaler.
# Tips      :
# Argument  : ref. of an array of numbers.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub sum_array{		# usage:  $output = &sum_array(\@any_array);
	my($sum);
	foreach $item(@{$_[0]}){
		$sum += $item;
	}
	\$sum;
}

#________________________________________________________________________
# Title     : print_seq_in_block
# Usage     : &print_seq_in_block (\$block_leng, 'i',\%h1, 'sort', \%h2, \%hash3,,,);
# Function  : gets many refs  for one scalar  or hashes and prints
#               the contents in lines of \$block_leng(the only scalar ref. given) char.
# Example   : If there are 3 hashes; (in the order of \%hash3, \%hash2, \%hash1)
#             >> 1st Hash        >> 2nd Hash         >> 3rd Hash          Output >>>
#             Name1  THIS-IS-    Name123  eHHHHHHH   Name123  12222223    Name123  12222223
#                                                              Name123  eHHHHHHH
#             options should be single chars                              Name1    THIS-IS-
# Warning   :
# Class     :
# Keywords  :
# Options   : 'o' or 'O' => ordered hash print, 'n' or'N' => no space between blocks.
#             's' or 'S' => printout sorted by seq names.(all options can be like \$sort
#             while $sort has 's' as value. naked number like 100 will be the
#             block_length. 'i' or 'I' => interlaced print.(this requires
#             identical names in hashes)
#             Example of ( no option, DEFAULT )        # Example of ( 'i' or 'I' option, INTERLACE )
#             6taa           ------ATPADWRSQSIY      #   6taa           ------ATPADWRSQSIY
#             2aaa           ------LSAASWRTQSIY      #   6taa           ------CCHHHHCCCCEE
#             1cdg           APDTSVSNKQNFSTDVIY      #   6taa           ------563640130000
#                                           #
#             6taa           ------CCHHHHCCCCEE      #   2aaa           ------LSAASWRTQSIY
#             2aaa           ------CCHHHHCCCCEE      #   2aaa           ------CCHHHHCCCCEE
#             1cdg           CCCCCCCCCCCCCCCCEE      #   2aaa           ------271760131000
#                                           #
#             6taa           ------563640130000      #   1cdg           APDTSVSNKQNFSTDVIY
#             2aaa           ------271760131000      #   1cdg           CCCCCCCCCCCCCCCCEE
#             1cdg           675232723600000000      #   1cdg           675232723600000000
#                                           #
#             Example of ('s' or 'S' option, SORT )    # Example of ('o' or 'O' option, ORDERED by input hashes
#             1cdg           APDTSVSNKQNFSTDVIY      #   6taa           ------ATPADWRSQSIY
#             2aaa           ------LSAASWRTQSIY      #   2aaa           ------LSAASWRTQSIY
#             6taa           ------ATPADWRSQSIY      #   1cdg           APDTSVSNKQNFSTDVIY
#                                           #
#             1cdg           CCCCCCCCCCCCCCCCEE      #   6taa           ------CCHHHHCCCCEE
#             2aaa           ------CCHHHHCCCCEE      #   2aaa           ------CCHHHHCCCCEE
#             6taa           ------CCHHHHCCCCEE      #   1cdg           CCCCCCCCCCCCCCCCEE
#                                           #
#             1cdg           675232723600000000      #   6taa           ------563640130000
#             2aaa           ------271760131000      #   2aaa           ------271760131000
#             6taa           ------563640130000      #   1cdg           675232723600000000
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : many refs  for hash (one for bottm, one for top, etc,top hash is usually
#               to denote certain caculations or results of the bottom one
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub print_seq_in_block{ my($k,$c,$i,$s,$t,@in,$intl,$z,$diff,$offset); my($bl)=60;
  $sort =0; $n_space =0; $ordered =0; $gap_char='-'; my($na,$larg,$names,$seq); my($n)=13;
  undef(@in);  sub numerically{ $a <=> $b; } my(@in_ar, $bl_passed);

  for($k=0; $k< @_ ;$k++){
	 if( !ref($_[$k]) ){    # when inputs are not ref.
		if($_[$k]=~ /^[\-]*([\d]{1,4})$/){ $bl=$1 if $1>0; $bl_passed=1 if $1>0; next;}  #<--   option handling
		if($_[$k]=~ /^[\-sS]$/){ $sort = 1; next;}
		if($_[$k]=~ /^[\-nN]$/){ $n_space = 1;next;}
		if($_[$k]=~ /^[\-iI]$/){ $intl = 1; next;}
		elsif($_[$k]=~ /^[\-oO]$/){ $ordered = 1;}
	 }
	 elsif( ref($_[$k]) eq "SCALAR" ){     #<--   option handling
		if(${$_[$k]}=~ /^[\-]*([\d]{1,4})$/){$bl=$1 if $1>0; $bl_passed=1 if $1>0; next;}      # the scalar input
		if(${$_[$k]}=~ /^[\-sS]$/){$sort = 1;next;}                 # should shorter than 5
		if(${$_[$k]}=~ /^[\-nN]$/){$n_space = 1;next;}              # to be recognized as
		if(${$_[$k]}=~ /^[\-iI]$/){$intl = 1;next;}                 # options, be it number or
		elsif(${$_[$k]}=~ /^[oO]$/){$ordered = 1;}                  # or chars.
	 }
	 elsif(ref($_[$k]) eq "HASH") { push(@in,  $_[$k]); } #<-- seqn handling hash
	 elsif(ref($_[$k]) eq "ARRAY"){ push(@in, &convert_arr_and_str_2_hash($_[$k], $k));} #<-- conv array to hash.
	 elsif(ref($_[$k]) eq "SCALAR"){ push(@in,&convert_arr_and_str_2_hash($_[$k], $k));} #<-- conv array to hash.
  }
  #########  HASH input handling ############
  for($k=0; $k<@in; $k++){
	   if(ref($in[$k]) eq "HASH"){ %{"input$k"}=%{$in[$k]};
	   if($sort == 1){
		  $keys_long=join("", keys(%{"input$k"}) );
	  if ($keys_long =~ /[\d\.]+/){
	     @{"names$k"}= sort numerically keys(%{"input$k"}); }   # numerical sort
	  elsif($keys_long =~ /[\w\.\,]+/){
	     @{"names$k"}= sort keys(%{"input$k"}); 	  } }       # normal sort
	  else{ @{"names$k"}= keys(%{"input$k"}); }
	   for($i=0; $i< @{"names$k"}; $i++){
		  if(${"input$k"}{${"names$k"}[$i]} =~ /\,/){               # remove ','
	     ${"input$k"}{${"names$k"}[$i]}=~ s/\,//g;}
	  }
	   }
  }

  for($z=0; $z < @in; $z++){
	for($t=0;$t<@{"names$z"};$t++){ $na=${"names$z"}[$t];$s=${"input$z"}{$na};
	   $larg=length($s) if length($s)> $larg;
	   $n=length($na) if length($na) > $n;
	   if($s =~ /\-/){ $gap_char='-'; }elsif( $s =~ /\./){  $gap_char='.';  }
	   if (length($s)<$larg){$offset=length($s);$diff=$larg-$offset; substr($s,$offset,$larg)="$gap_char"x$diff;}}}
  if($ordered== 1){  $bl=$larg if (($larg < 60)&&($bl_passed != 1));
	  for($c=0; $c < @in; $c++){
		 for($k=0; $k < $larg; $k+=$bl){
	    for($i=0; $i <=$#{"names$c"}; $i++){
	      $names= ${"names$c"}[$i];
	      $seq= substr( ${"input$c"}{$names},  $k,  $bl);
	      printf ("%-$n s  %-$bl s\n", $names,$seq);
	    }
	    print "\n" unless($n_space ==1);
	 }print "\n";
	  }
  }
  elsif($intl==1){   $bl=$larg  if (($larg < 50)&&($bl_passed != 1));
		for($k=0; $k < $larg; $k+=$bl){
		  for($i=0; $i < @{"names0"}; $i++){
			for($c=0; $c <= $#in; $c++){ $names=${"names$c"}[$i]; $seq=substr(${"input$c"}{$names}, $k, $bl);
	       printf ("%-$n s  %-$bl s\n", $names,$seq); } print "\n" unless($n_space ==1);} print "\n";} print "\n" unless($n_space ==1);}
  else{ for($k=0; $k < $larg; $k+=$bl){   $bl=$larg if (($larg < 50)&&($bl_passed != 1));
		  for($c=0; $c < @in; $c++){  # $n is the name space size
	    for($i=0; $i < @{"names$c"}; $i++){ $names=${"names$c"}[$i]; $seq=substr(${"input$c"}{$names},$k,$bl);
	    printf ("%-$n s  %-$bl s\n", $names,$seq);} print "\n" unless($n_space ==1);}print "\n";}}
}

#________________________________________________________________________
# Title     : cls
# Usage     : &cls;
# Function  : clears screen
# Example   : &cls;
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub cls{   my($cls) = `clear`;  print $cls;  }


#________________________________________________________________________
# Title     : read_dir_names_only
# Usage     : @all_dirs_list = @{&read_dir_names_only(\$absolute_path_dir_name, ....)};
# Function  : read any dir names and and then put in array.
# Example   :
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Todo      :
# Author    : A Biomatic
# Version   : 3.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_dir_names_only{
  my($in_dir, $i,$k, @possible_dirs,
	  @final_files, $full_dir, $pwd, $path,@read_files);
  $pwd=`pwd`; chomp($pwd); $full_dir=1;
  for($k=0; $k < @_; $k++){
	 if   ( ($_[$k] eq '.') || !(defined($_[$k]))){  $in_dir=$pwd;  }
	 elsif(!(ref($_[$k]))){   $in_dir=$_[$k];   }
	 elsif(ref($_[$k])){      $in_dir =${$_[$k]};    }
	 if($in_dir =~ /^([\w\-\.]+)$/){  $in_dir="$pwd\/$in_dir"; $full_dir = 0; }
	 else{ $full_dir =1; }
	 ##########  Main READING PART ##########
	 opendir(DIR1,"$in_dir");
	 @read_files = readdir(DIR1);
	 for($i=0; $i < @read_files; $i ++){
		$read_files[$i]="$in_dir\/$read_files[$i]";
		if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
		  $read_files[$i]=~s/\.\///; ## removing ./ in front of dirs (in bash)
		  push(@final_files, "$read_files[$i]");
		}
	 }
  }
  return([sort @final_files]);
}
