/*

実行の際は以下を打ち込む
gcc -c spectrogram.c -lm
gcc fft.c -lm -c
gcc spectrogram.o fft.o -lm -o spectrogram
./spectrogram a00.raw out.dat

gnuplotのコマンド
set pm3d map
splot "hogehoge.dat"

*/
/*

注意

音ファイルのbit数に合わせて、short型、double型を変更するのを忘れずに
*/

#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>

#define L 320 //標本数( t[ms]/(1/x[kHz]) )
#define　sample_L 0.02 //標本間隔(t秒)実験では0.02秒間隔だった
#define OverRap 0.01 //オーバーラップする時間(t秒) 実験では0.01秒間隔だった
#define N 1024 //データ点数
#define frequency_interval 0.015625 //出力の際の周波数間隔(x[kHz]/1024[FFT長])
#define bit 16 //ビット数


double Euclid(double r[], double unknown[]);
void powerspectrum(char *arg, double S[]);
int getPowerSpectrum(char *arg, double S[], int position);
//int getPowerSpectrum(short y[], char *arg, double S[], int position);

int main(int argc, char *argv[])
{
   /* argv[0]はプログラム名 */
   /* argv[1]は第1引数名で未知の音声データ */
   /* argv[2]は第2引数名で出力データ */

  FILE *fp;

  int n,m,i;
  //最小ユークリッド距離の番号を保存
  int min=1;
  //認識判定結果を保持する変数
  int run=0;
  //フレームの数
 int frame = 0;

  //未知のデータ
  double unknown[N] = {0};
//周波数軸
  double fre=0;
  //時間
  double time=0;

  //stat関数を使用する際にファイルの状態を保存するために使われる構造体
  struct stat filesize;
  //標本数を保持(バイト / 2)
  int sample;
  int position = 0;
  int byte = 0;
  
  //バイトの計算
  byte = bit/8;
  
  //バイナリ形式で読み込み
  if((fp = fopen( argv[2], "w" ))==NULL)
    {
      printf("ファイルを開くことができません\n");
      exit(1);
    }
  

  //ファイルの状態をstatによって構造体に保存
  //http://www.c-lang.net/stat/index.html
  stat(argv[1], &filesize);
  //「.st_size」を使い総標本数を求める
  sample = filesize.st_size / byte;
  //読み込む長さを計算
  frame = (int) sample / (L/2);
  //総バイト数分配列を保持、票本数以上配列を取得
  //y = (short *)malloc(filesize.st_size);



      
  //標本の読み込む位置を計算しつつユークリッド距離を計算
  for(m=0; m<frame; m++)
    {     
     //フレームの位置を計算、バイトに直すため2倍
      position = (m*L/2)*2;
      //run = getPowerSpectrum(y, argv[6], unknown, position);
      run = getPowerSpectrum(argv[1], unknown, position);
      

	  for(i=0; i<N/2; i++){
	      fprintf(fp, "%f\t%f\t%f\n", time, fre, unknown[i]);
	      //周波数間隔は標本数によって変わるので注意
	      fre += frequency_interval;
	  }
	  
	  fprintf(fp, "\n");
	  fre = 0;
	  time += OverRap;
    }

fclose(fp);

  return 0;
}

//int getPowerSpectrum(short y[], char *arg, double S[], int position){
int getPowerSpectrum(char *arg, double S[], int position){

  //stat関数を使用する際にファイルの状態を保存するために使われる構造体
  struct stat filesize;

  FILE *fp;
  int i;

  //総標本数を保存する変数
  int sample;
  //一度ダブル型で計算するために用意された変数
  double Hamming = 0;
  //フーリエ変換用の配列
  double yr[N] = {0};
  double yi[N] = {0};
  //母音の認識をするかどうかの判別に用いる変数
  double Decision=0;
  int byte2 = 0;
  
  //バイトの計算
  byte2 = bit/8;
  
    //バイト数分、動的メモリ確保する変数
  short *y;

  //バイナリ形式で読み込み
  if((fp = fopen( arg, "rb" ))==NULL)
    {
      printf("ファイルを開くことができません\n");
      exit(1);
    }
    
    
  //ファイルの状態をstatによって構造体に保存
  //http://www.c-lang.net/stat/index.html
  stat(arg, &filesize);
  //「.st_size」を使い総標本数を求める
  sample = filesize.st_size / byte2;
  //総バイト数分配列を保持、票本数以上配列を取得
  y = (short *)malloc(filesize.st_size);


  //fseekより指定された位置から読み込むようにセット
  fseek(fp, position, SEEK_SET);
  fread( y, byte2, L, fp );

  Decision=0;
  for(i=0;i<L;i++)
  {
     yi[i]=0;
     yr[i]=0;
  }
  
  //窓掛け
  for(i=0; i<L; i++)
    {
      //判定のためにここで計算
      Decision += y[i]*y[i];

      Hamming = ( 0.5 - 0.5 * cos( (2*i*M_PI) / (L - 1) ) ) * y[i];
      yr[i] = Hamming;
    }

  Hamming=0;
  free(y);

  //判定
  //しきい値は300に設定するともっとも結果が正しいものとなった。
  //Decision = sqrt( Decision / L );
  //フーリエ変換
      fft(N, yr, yi); 
      
      //パワースペクトル計算
      for(i=0; i<N; i++)
	{
	  S[i] = 10 * log10( yi[i]*yi[i] + yr[i]*yr[i] ) +60; 
	  
	  //printf("%f\t%f\n", fre, S[i]);
	  //fre += 0.015625;
	}

      fclose(fp);
      return 1;
}
