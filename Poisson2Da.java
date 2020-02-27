package poisson2da;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

// Program 2D Poisson equation 


public class Poisson2Da {
	public static void main(String... args) {
		int N = 100;
		// Constants
		double X = 1.0; //計算領域の大きさ
		double e0 = 8.85e-12;//真空の誘電率
		int center = N/2;//中心の座標
		double delta = X/N; //格子の間隔
		double Conv = 1.0e-6;//収束と半手する前回ループとの差
		
		double[][] phi = new double[N][N]; // 計算すべき電位
		double[][] rho = new double[N][N]; // 電荷密度
		double MaxPhi;//最大電位
		double MaxErr;//最大のエラー
		double CurErr;//現在のエラー
		double Prev_phi;//前のループのphi
		double Ex,Ey;//電場
//		int i,j;
		int loop;//繰り返しのカウンタ
//		PrintWriter fileWriter=null; //出力ファイルライター
		
		// initialize phi,rho
		for(int i = 0;i<N;i++) {
			for(int j = 0;j<N;j++) {
				phi[i][j] = 0.0;
				rho[i][j] = 0.0;
			
			}
		}
		
		// 電荷を置く 半径５ｃｍ、電荷密度 1x10^-8[C/m2]
		for(int i = 0;i<N;i++) {
			for(int j = 0;j<N;j++) {
				if(((center-i)*(center-i) + (center-j)*(center-j))*delta*delta<0.05*0.05) {
					rho[i][j] = 1.0e-8;
				}
			}
		}
		
		// 繰り返し計算
		loop = 0;
		MaxPhi = 1.0e-10; //系内の最大の電位を入れる変数。ある有限の値を入れておく（ゼロ割防止）
		
		do {
			if(loop%1000 == 0) {
				System.out.printf("%05d  %e\n", loop,MaxPhi);
			}
			MaxErr = 0.;
			CurErr = 0.;
			for(int i=1;i<N-1;i++) {
				for(int j=1;j<N-1;j++) {
					Prev_phi = phi[i][j]; // 前回のループのphiを入れておく
					phi[i][j] = 0.25 * (rho[i][j]*delta*delta/e0 + phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]);
										//Poisson方程式のphiを計算
					if(MaxPhi<Math.abs(phi[i][j])) {
						MaxPhi = phi[i][j]; // 最大の電位が更新されたら書き換え
					}
					CurErr= (Math.abs(phi[i][j] - Prev_phi))/MaxPhi; //前回ループと新しい答えの差をMaxPhiで規格化
					if(MaxErr<CurErr) {
						MaxErr = CurErr; // 誤差の最大を常にMaxErrに持つようにする。
					}
				}
				
			}
			loop++;
		}while(MaxErr>Conv);
		
		//ポテンシャルの出力
		PrintWriter fileWriter=null; //出力ファイルライター
		File file = new File("Outputfile/phi.avd");
		try {
			fileWriter = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			for(int i = 0;i<N;i++) {
				for(int j = 0;j<N;j++) {
			fileWriter.printf("%e %e %e\n", delta*i,delta*j,phi[i][j]);
				}
			}
		}catch (IOException e) {
			e.printStackTrace();
		}finally {
			if(fileWriter != null) {
				try {
					fileWriter.close();
				}catch (Exception e) {
					// TODO: handle exception
				}
			}
		}
		
		//電場の出力
	}
	
}
