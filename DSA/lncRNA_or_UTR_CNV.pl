
$start_time = time();
#测试修改代码
open(INPUT,"D:/sSNP/Se_Data/SNP_TCGAsample.txt")||die "can't open:$!\n";
while(defined($aaa=<INPUT>)){
   chomp($aaa);
   $aaa=substr($aaa,0,16);
   $SM_hash{$aaa}=1;
}
close(INPUT);
open(INPUT,"D:/sSNP/Data/TCGA_33_multiple_cancer/annotation/gdc_sample_sheet_all.tsv")||die "can't open:$!\n";
$aaa=<INPUT>;
while(defined($aaa=<INPUT>)){
   chomp($aaa);
   @arr=split(/\t|\.grch38/,$aaa);
   $folder_hash{$arr[1]}=$arr[7];
}
close(INPUT);

my ($dirname, $endir);
$dirname = "D:/sSNP/Data/TCGA_33_multiple_cancer/GISTIC/";
opendir (DIR, $dirname ) || die "Error in opening dir $endir\n";
while( ($endir = readdir(DIR))){
     if (!($endir=~/\.|_/)){

open(INPUT,"D:/sSNP/Data/TCGA_33_multiple_cancer/GISTIC/".$endir."/all_lesions.conf_95.txt")||die "can't open:$!\n";
$aaa=<INPUT>;
chomp($aaa);
@arr=split(/\t/,$aaa);
my @sample_col;
my $n=0;
my $first_line="gene";
foreach $i(@arr){
   if($folder_hash{$i}){
      if($SM_hash{$folder_hash{$i}}){
         $first_line=$first_line."\t$folder_hash{$i}";
         push (@sample_col,$n);
         delete $SM_hash{$folder_hash{$i}};
      }
   }
   $n++;
}
if(@sample_col){
 while(defined($aaa=<INPUT>)){
   chomp($aaa);
   @arr=split(/\t/,$aaa);
   if($arr[0]=~/CN/){
      @posi=split(/\(|:|-/,$arr[2]);
      if(!$hash{$posi[0]."_".$posi[1]."_".$posi[2]}){
        foreach $i(@sample_col){
          if($hash{$posi[0]."_".$posi[1]."_".$posi[2]}){
             $hash{$posi[0]."_".$posi[1]."_".$posi[2]}=$hash{$posi[0]."_".$posi[1]."_".$posi[2]}."_".$arr[$i];
          }
          else{
             $hash{$posi[0]."_".$posi[1]."_".$posi[2]}=$arr[$i];
          }
        }
      }
   }
 }
}
close(INPUT);

if(@sample_col){
mkdir("D:/sSNP/Data/TCGA_23_cancer/UTR_cnv/$endir");
open(OUTPUT,">D:/sSNP/Data/TCGA_23_cancer/UTR_cnv/$endir/UTR_cnv.txt")||die "can't open:$!\n";
print OUTPUT  $first_line."\n";
open(INPUT,"D:/sSNP/Se_Data/gencode_UTR_genePosition.txt")||die "can't open:$!\n";
while(defined($aaa=<INPUT>)){
   chomp($aaa);
   @arr=split(/\t/,$aaa);
   my $UTR_cover=0;
   my $over_size=0;
   my $UTR_cnv;
   foreach $key(keys %hash){
      @Region_posi=split(/_/,$key);
      if($arr[0] eq $Region_posi[0]){
          $over_size=Overlop_size($arr[1],$arr[2],$key);
          if($over_size>$UTR_cover){
              $UTR_cover=$over_size;
              $UTR_cnv=$hash{$key};
          }
      }
   }
   if($UTR_cover>0){
      print OUTPUT $arr[-1]."\t";
      @out=split(/_/,$UTR_cnv);
      foreach $i(@out){
         print OUTPUT $i."\t";
      }
      print OUTPUT "\n";
   }
}
close(INPUT);
close(OUTPUT);
}

if(@sample_col){
mkdir("D:/sSNP/Data/TCGA_23_cancer/lncRNA_cnv/$endir");
open(OUTPUT,">D:/sSNP/Data/TCGA_23_cancer/lncRNA_cnv/$endir/lncRNA_cnv.txt")||die "can't open:$!\n";
print OUTPUT  $first_line."\n";
open(INPUT,"D:/sSNP/Se_Data/lncRNA_genePosition.txt")||die "can't open:$!\n";
while(defined($aaa=<INPUT>)){
   chomp($aaa);
   @arr=split(/\s/,$aaa);
   my $lncRNA_cover=0;
   my $over_size=0;
   my $lncRNA_cnv;
   foreach $key(keys %hash){
      @Region_posi=split(/_/,$key);
      if($arr[0] eq $Region_posi[0]){
              my $single_lnc_size=0;
              for(my $i=1;$i<$#arr;$i=($i+3)){
            $single_lnc_size=Overlop_size($arr[$i],$arr[$i+1],$key)+$single_lnc_size;
                  }
                  $over_size=$single_lnc_size;
          if($over_size>$lncRNA_cover){
              $lncRNA_cover=$over_size;
              $lncRNA_cnv=$hash{$key};
          }
      }
   }
   if($lncRNA_cover>0){
      print OUTPUT $arr[-1]."\t";
      @out=split(/_/,$lncRNA_cnv);
      foreach $i(@out){
         print OUTPUT $i."\t";
      }
      print OUTPUT "\n";
   }
}
close(INPUT);
close(OUTPUT);
}

}
}
closedir(DIR);

sub Overlop_size{
   my $size=0;
   @v=split(/_/,$_[2]);
   if($_[0]>$v[2]){}
   if($_[1]<$v[1]){}
   if($_[0]<$v[1]){
      if($_[1]>=$v[1]){
         if($_[1]<=$v[2]){
            $size=$_[1]-$v[1];
         }
      }
   }
   if($_[0]>=$v[1]){
      if($_[0]<=$v[2]){
         if($_[1]>$v[2]){
            $size=$v[2]-$_[0];
         }
      }
   }
   if($_[0]>=$v[1]){
      if($_[1]<=$v[2]){
         $size=$_[1]-$_[0];
      }
   }
   if($_[0]<=$v[1]){
      if($_[1]>=$v[2]){
         $size=$v[2]-$v[1];
      }
   }
   return $size;
}

$end_time = time();
$elapsed_time = $end_time - $start_time;
print $elapsed_time;