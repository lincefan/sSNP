#!/usr/bin/perl
#111111
$start_time = time();

open(INPUT,"D:/sSNP/priUTR_position/Somatic_mutations_TCGA_930_mt3sample.txt")||die "can't open:$!\n";
while(defined($aaa=<INPUT>)){
     chomp($aaa);
     @arr=split(/\t/,$aaa);
     $hash{$arr[0]}=$arr[0];
     push @{$hash{$arr[0]}},$arr[1] ;##take attention to the form @{$hash{$key}} and ${$hash{$key}};
}
close(INPUT);
foreach $key(keys %hash){
     @{$hash{$key}}=sort {$a<=>$b}(@{$hash{$key}}); ##<=> means sorting by numeric order
}

open(OUTPUT,">D:/sSNP/priUTR_position/priUTR_poisi_mt3.txt")||die "can't open:$!\n";
$dirname = "D:/sSNP/priUTR_position/";
opendir (DIR, $dirname ) || die "Error in opening dir $dirname\n";
while( ($filename = readdir(DIR))){
if ($filename=~"chr"){

 open(INPUT_1,"D:/sSNP/priUTR_position/".$filename)||die "can't open:$!\n";
 while(defined($bbb=<INPUT_1>)){
        chomp($bbb);
        @brr=split(/\t/,$bbb);
        @me=split(/\;/,$brr[4]);
        shift(@me);
        my $siz=@me;
        my $chara;
        my $keep=0;
        my $inte_1;
        for($i=0;$i<$siz;$i=$i+2){
           @inte=split(/_/,$me[$i+1]);
           $inte_1=find_SM($inte[0],$inte[1],@{$hash{$me[$i]}});
           $chara=$chara.";".$inte_1;
           if($inte_1){
              $keep=1;
           }
        }
        if($chara ne ";"){
         print OUTPUT $bbb."\t".$chara."\n";
        }

 }
close(INPUT_1);

}
}
closedir (DIR);
close(OUTPUT);

$end_time = time();
$elapsed_time = $end_time - $start_time;
print $elapsed_time;

sub find_SM{
   my $o_SM;
   my $in_start=@_[0];
   my $in_end=@_[1];
   my @arr=@_[2..(@_-1)];
   if($in_start<$arr[$#arr] && $in_end>$arr[0]){
      $index_start=bfind_start($in_start,0,$#arr,@arr);
      $index_end=bfind_end($in_end,0,$#arr,@arr);
      foreach $it(@arr[$index_start..$index_end]){
         $o_SM=$o_SM.";".$it;
      }
   }
   return $o_SM;
}
sub bfind_start{  #shouldn't change order  when implement pe($s,@arr)
   my $n=@_[0]; #my pramater only use in subway
   my $start=@_[1];
   my $end=@_[2];
   my @arr=@_[3..(@_-1)]; #!the second pramater use as array
   my $mid=int(($end+$start)/2);
   if($n==$arr[$start]){  return $start;}
   elsif($n==$arr[$end]){ return $end;}
   elsif(($end-$start)==1){ return $end;}
   elsif($n==$arr[$mid]){ return $mid;}
   elsif($n>$arr[$mid]){
      $start=$mid;
      bfind_start($n,$start,$end,@arr);
   }
   elsif($n<$arr[$mid]){
      $end=$mid;
      bfind_start($n,$start,$end,@arr);
   }
}
sub bfind_end{  #shouldn't change order  when implement pe($s,@arr)
   my $n=@_[0]; #my pramater only use in subway
   my $start=@_[1];
   my $end=@_[2];
   my @arr=@_[3..(@_-1)]; #!the second pramater use as array
   my $mid=int(($end+$start+1)/2);
   if($n==$arr[$start]){  return $start;}
   elsif($n==$arr[$end]){ return $end;}
   elsif(($end-$start)==1){ return $start;}
   elsif($n==$arr[$mid]){ return $mid;}
   elsif($n>$arr[$mid]){
      $start=$mid;
      bfind_end($n,$start,$end,@arr);
   }
   elsif($n<$arr[$mid]){
      $end=$mid;
      bfind_end($n,$start,$end,@arr);
   }
}