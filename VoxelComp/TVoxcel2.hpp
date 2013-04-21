//ボクセル系の関数
//ボクセル構造体の要素に、「内部に含まれる点の個数」「信頼度= (ボクセル内部の個数 / 全体の個数) * (８＾階層レベル)」

/*mask読み込み*/
// #include <windows.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//ボクセルの分割数(8分木構造体なら8で固定)
#define CHILD_VOXCEL_NUM 8


/*----------------------------------------------------------------------------
	ボクセル構造体
----------------------------------------------------------------------------*/
struct TVoxcel {
  public:
    int id;						//識別子
    TVoxcel *parent;			//親ノード（これがNullのとき、このノードは根）
    TVoxcel **child;			//子ノード（これがNullのとき、このノードは葉）
    int child_size;				//子ノードの数:8分木なら基本8
    int level;					//階層レベル（これが０のとき、このノードは根）
    double min_x, min_y, min_z;	//この領域の最小の点
    double max_x, max_y, max_z;	//この領域の最大の点
    int point_number;			//内包する点群の個数
    double value;				//このボクセルが持つ値(=階層数)
    double state;				//このボクセルの状態(-1:存在しない,0:処理中,1:存在,2:固定、など)
};

TVoxcel *root;						//八分木作成用ルートボクセル

//プロトタイプ宣言
void Rec_deleteTVoxcel(TVoxcel *voxcel);


/*----------------------------------------------------------------------------
	ボクセル初期化
	--------------------------------------------------------------------------
	child_sizeが8の場合のみ8分木構造と判断し、
	min_x,min_y,min_z,max_x,max_y,max_zにparentから適切な値を設定する。
	それ以外の場合は初期値として0を代入して返す。
----------------------------------------------------------------------------*/
//voxcel: 初期化するボクセル
//parent: 親ボクセル
//child_num: parentに対するvoxcelの番号(何番目の子か)
//state: ボクセルの状態
void initTVoxcel(TVoxcel *voxcel, TVoxcel *parent = 0, int child_num = 0, double state = -1) {
	if (parent != 0) {
		voxcel->level = parent->level + 1;
		voxcel->child_size = parent->child_size;  //親と同じ分木が入る：基本は８
		voxcel->id = static_cast<int>((parent->id + pow(1.0 * voxcel->child_size, 1.0 * voxcel->level) * (child_num + 1)));  //疑似8進法:000000000
		voxcel->parent= parent;
		voxcel->child =NULL;
		if (voxcel->child_size == 8) {  //8分木時：条件演算子→代入される物　＝　条件式　？　真：偽
			voxcel->min_x = (child_num % 2 == 0)     ? parent->min_x : (parent->min_x + parent->max_x) / 2;
			voxcel->min_y = (child_num / 2 % 2 == 0) ? parent->min_y : (parent->min_y + parent->max_y) / 2;
			voxcel->min_z = (child_num / 4 == 0)     ? parent->min_z : (parent->min_z + parent->max_z) / 2;
			voxcel->max_x = (child_num % 2 == 0)     ? (parent->min_x + parent->max_x) / 2 : parent->max_x;
			voxcel->max_y = (child_num / 2 % 2 == 0) ? (parent->min_y + parent->max_y) / 2 : parent->max_y;
			voxcel->max_z = (child_num / 4 == 0)     ? (parent->min_z + parent->max_z) / 2 : parent->max_z;
		}
    else{  //八分木ではない時
			voxcel->min_x = 0, voxcel->min_y = 0, voxcel->min_z = 0;
			voxcel->max_x = 1, voxcel->max_y = 1, voxcel->max_z = 1;
		}
		voxcel->point_number=0;
		voxcel->value = 0;
		voxcel->state = state;
  }
  else {  //ルートボクセル用
		voxcel->level = 0;
		voxcel->id    = 0;
		voxcel->parent= NULL;
		voxcel->child = NULL;
		voxcel->child_size = CHILD_VOXCEL_NUM;
		voxcel->min_x = 0;
		voxcel->min_y = 0;
		voxcel->min_z = 0;
		voxcel->max_x = 1;
		voxcel->max_y = 1;
		voxcel->max_z = 1;
		voxcel->point_number = 0;
		voxcel->value = 0;
		voxcel->state = state;
	}
}

/*----------------------------------------------------------------------------
	ボクセルを分割
	--------------------------------------------------------------------------
	ボクセルの分割処理を行い，生成した子要素ボクセルを初期化する
----------------------------------------------------------------------------*/
//voxcel: 分割するボクセル
//return: 分割処理を行ったか？(voxcelが既に分割されている場合はfalseを返す)
bool devideTVoxcel(TVoxcel *voxcel) {
	if (voxcel->child != 0) {
    return false;
  }
	voxcel->child = new TVoxcel * [voxcel->child_size];  //まだ枝があるということで、8個作っておく。
	for (int i = 0; i < voxcel->child_size; i++){  //１〜８(8分木なら8、4分木なら4が入る)
		voxcel->child[i] = new TVoxcel();
		initTVoxcel(voxcel->child[i], voxcel, i);
	}
	return true;
}

/*----------------------------------------------------------------------------
	木構造ボクセルの削除
	--------------------------------------------------------------------------
	ボクセルを親から削除し、自分の子要素もすべて消す
----------------------------------------------------------------------------*/
//voxcel: 削除するボクセル
void deleteTVoxcel(TVoxcel *voxcel) {
	if(voxcel->child != 0) {
		fprintf(stderr, "ボクセル削除[deleteTVoxcel]");  //確認用
    Rec_deleteTVoxcel(voxcel);
	}
}


/*----------------------------------------------------------------------------
	指定した座標をボクセルが含むか？
	--------------------------------------------------------------------------
	引数(x,y,z)で指定した領域をボクセルが含むか判定
----------------------------------------------------------------------------*/
//voxcel: ボクセル
//x,y,z:  3次元座標値
//return: 含むならtrue
bool includePointInTVoxcel(TVoxcel *voxcel, double x, double y, double z) {
	if (voxcel->min_x <= x
      && voxcel->max_x >= x
      && voxcel->min_y <= y
      && voxcel->max_y >= y
      && voxcel->min_z <= z
      && voxcel->max_z >= z) {
    voxcel->point_number++;
    return true;
	}
  else {
    return false;
  }
}

/*----------------------------------------------------------------------------
	指定した座標を含む指定階層のボクセルを返す
	--------------------------------------------------------------------------
	引数(x,y,z)で指定した領域を含むボクセルを返す
	返すボクセルの階層数はlevelで指定し、
	値がマイナスの場合は子要素を持たないボクセルを返す
----------------------------------------------------------------------------*/
//voxcel: ボクセル
//x,y,z:  3次元座標値
//return: 含むならtrue
TVoxcel *getTVoxcel(TVoxcel *voxcel, double x, double y, double z, int level = -1) {
	TVoxcel *tmp_vox = voxcel;	// 現在返す候補のボクセル
	while (tmp_vox->child != 0 && tmp_vox->level != level) {
		for (int i = 0; i < tmp_vox->child_size; i++) {
			if (tmp_vox->child[i]->min_x <= x
          && tmp_vox->child[i]->max_x >= x
          && tmp_vox->child[i]->min_y <= y
          && tmp_vox->child[i]->max_y >= y
          && tmp_vox->child[i]->min_z <= z
          && tmp_vox->child[i]->max_z >= z) {
        tmp_vox = tmp_vox->child[i];
				break;
			}
			if (i == tmp_vox->child_size - 1) {
        return 0;
      }
		}
	}
	if (tmp_vox->level == level || level < 0) {
    return tmp_vox;
  }
	else {
    return 0;
  }
}

/////////////////////////////////////////////////////////////////////////////////
//以下は汎用ではない
/////////////////////////////////////////////////////////////////////////////////


/*----------------------------------------------------------------------------
	ボクセルを指定階層まで分割：cppファイルから呼ばれる
	--------------------------------------------------------------------------
	引数(x,y,z)で指定した領域を含むボクセルを指定階層まで分割する
----------------------------------------------------------------------------*/
//voxcel: ボクセル
//x,y,z:  3次元座標値
//level:  指定階層
//state:  分割してできたボクセルの状態(前回の状態より大きな状態であれば値は更新される)
void devideTVoxcelByPoint(TVoxcel *voxcel, double x, double y, double z, int level, double state = 1.0) {
	if (includePointInTVoxcel(voxcel, x,y,z)) {
		if (level>voxcel->level) {
			if (voxcel->child == 0) {
        devideTVoxcel(voxcel);
      }
			for (int i = 0; i < voxcel->child_size; i++) {
        devideTVoxcelByPoint(voxcel->child[i], x, y, z, level, state);
      }
		}
		voxcel->value = voxcel->level;
		if (state > voxcel->state) {
      voxcel->state = state;
    }
	}
}


/*----------------------------------------------------------------------------
	ボクセルの状態を変更
	--------------------------------------------------------------------------
	引数pre_stateの状態にあるボクセルすべての値を,set_stateに置き換える
----------------------------------------------------------------------------*/
//voxcel		ボクセル
//pre_state		置き換えるボクセルの状態値
//set_state		変更後の状態の値
void replacedTVoxcelState(TVoxcel *voxcel, double pre_state, double set_state){
	if(voxcel->state==pre_state)voxcel->state=set_state;
	if(voxcel->child!=NULL)
		for(int i=0;i<voxcel->child_size;i++)
			replacedTVoxcelState(voxcel->child[i],pre_state,set_state);	
}

/*----------------------------------------------------------------------------
	ボクセルを削除(サブルーチン)
	--------------------------------------------------------------------------
	ボクセルを子要素を消し，自分も消す．
----------------------------------------------------------------------------*/
//voxcel: 削除するボクセル：cppからはdeleteを経由してrootが送られる
void Rec_deleteTVoxcel(TVoxcel *voxcel){
	if(voxcel->child != NULL){
		for(int i=0; i<8; i++){
			//TVoxcel* voxcel = voxcel->child[i];	
			//Rec_deleteTVoxcel(voxcel);  // 子以降を削除します
			if(voxcel->child[i] != NULL)Rec_deleteTVoxcel(voxcel->child[i]);  // 子以降を削除します
		}
	}
	if(voxcel->id < 80)	fprintf(stderr, "->");//確認用
	delete voxcel;  // 最後に自分を削除
}


/*----------------------------------------------------------------------------
	縮退処理において注目ボクセルを削除すると26近傍で穴があくかのチェック
	--------------------------------------------------------------------------
	引数filterの配列に対し、穴のあいた箇所を指定の条件ですべて通ることができるか
	通ることができるのは値が2の要素で、なおかつ注目ボクセルや端の8箇所のボクセル
	を除く箇所である。すべて通るならば、値が-1の箇所はすべて0に書き換えられる
----------------------------------------------------------------------------*/
//filter	フィルタ
//x,y,z		現在の位置
void erodeHoleCheck(int *filter, int x, int y, int z){
	if(filter[z*9+y*3+x]==-1) filter[z*9+y*3+x]=0;
	else if(filter[z*9+y*3+x]==2) filter[z*9+y*3+x]=3;
	else return;
	for(int k=-1;k<=1;k++)
		for(int j=-1;j<=1;j++)
			for(int i=-1;i<=1;i++)
				if(((z+k)*9+(y+j)*3+(x+i))>=0 && ((z+k)*9+(y+j)*3+(x+i))<27 &&!(
						((z+k)*9+(y+j)*3+(x+i))==0 ||
						((z+k)*9+(y+j)*3+(x+i))==2 ||
						((z+k)*9+(y+j)*3+(x+i))==6 ||
						((z+k)*9+(y+j)*3+(x+i))==8 ||
						((z+k)*9+(y+j)*3+(x+i))==18||
						((z+k)*9+(y+j)*3+(x+i))==20||
						((z+k)*9+(y+j)*3+(x+i))==24||
						((z+k)*9+(y+j)*3+(x+i))==26) )
							erodeHoleCheck(filter,x+i,y+j,z+k);
}

/*----------------------------------------------------------------------------
	縮退処理において注目ボクセルを削除しても26近傍で連結しているかのチェック
	--------------------------------------------------------------------------
	引数filterの配列に対し、値が2か3の要素を26近傍の移動ですべて通ることができるか
	すべて通るならば、値が2か3の箇所はすべて4に書き換えられる
----------------------------------------------------------------------------*/
//filter	フィルタ
//x,y,z		現在の位置
void erodeConnectionCheck(int *filter, int x, int y, int z){
	if(filter[z*9+y*3+x]==2 || filter[z*9+y*3+x]==3) filter[z*9+y*3+x]=4;
	else return;
	for(int k=-1;k<=1;k++)
		for(int j=-1;j<=1;j++)
			for(int i=-1;i<=1;i++)
				if(((z+k)*9+(y+j)*3+(x+i))>=0 && ((z+k)*9+(y+j)*3+(x+i))<27)
					erodeConnectionCheck(filter,x+i,y+j,z+k);
}

/*----------------------------------------------------------------------------
	ボクセルを指定階層において縮退処理
	--------------------------------------------------------------------------
	引数levelで指定した階層数において、縮退処理を実行する
	この処理で新しく生成されたボクセルの状態は0となる（通常可視のボクセルは1）
----------------------------------------------------------------------------*/
//voxcel	ボクセル
//level		指定階層
//root_vox	ルートボクセル
void erodeTVoxcel(TVoxcel *voxcel, int level, TVoxcel *root_vox=NULL){
	if(root_vox==NULL)root_vox=voxcel;
	if(voxcel->level==level){
		if(voxcel->state==1){
			int filter[27]; TVoxcel *tmp_vox; bool deleteFlag=false;
			double vox_pos[3]; double vox_length[3];
			vox_pos[0]=(voxcel->max_x+voxcel->min_x)/2;
			vox_pos[1]=(voxcel->max_y+voxcel->min_y)/2;
			vox_pos[2]=(voxcel->max_z+voxcel->min_z)/2;
			vox_length[0]=voxcel->max_x-voxcel->min_x;
			vox_length[1]=voxcel->max_y-voxcel->min_y;
			vox_length[2]=voxcel->max_z-voxcel->min_z;
			//printf("%f,%f,%f(%f,%f,%f)\n",vox_pos[0],vox_pos[1],vox_pos[2],vox_length[0],vox_length[1],vox_length[2]);
			//注目ボクセルの周囲のボクセル状況を記録
			//-1:6近傍のボクセルなし,0:26近傍のボクセルなしか注目ボクセル
			//1:ボクセル不可視,2:ボクセル存在
			for(int z=-1;z<=1;z++){
				for(int y=-1;y<=1;y++){
					for(int x=-1;x<=1;x++){
						tmp_vox = getTVoxcel(root_vox,
								vox_pos[0]+vox_length[0]*x,
								vox_pos[1]+vox_length[1]*y,
								vox_pos[2]+vox_length[2]*z,
								level);
						if(tmp_vox!=NULL && (9*(z+1)+3*(y+1)+(x+1))!=13 && tmp_vox->state!=-1){
							if(tmp_vox->state==0){
								filter[9*(z+1)+3*(y+1)+(x+1)]=1;
							}else{
								filter[9*(z+1)+3*(y+1)+(x+1)]=2;
							}
						}else{
							//printf("%d\n",(9*(z+1)+3*(y+1)+(x+1)));
							if( (9*(z+1)+3*(y+1)+(x+1))==4 ||
								(9*(z+1)+3*(y+1)+(x+1))==10||
								(9*(z+1)+3*(y+1)+(x+1))==12||
								(9*(z+1)+3*(y+1)+(x+1))==14||
								(9*(z+1)+3*(y+1)+(x+1))==16||
								(9*(z+1)+3*(y+1)+(x+1))==22){
									filter[9*(z+1)+3*(y+1)+(x+1)]=-1;
							}else{
								filter[9*(z+1)+3*(y+1)+(x+1)]=0;
							}
						}
					}
				}
			}
			//printf("%d,%d,%d %d,%d,%d %d,%d,%d\n",filter[9*0+3*0+0],filter[9*0+3*0+1],filter[9*0+3*0+2],
			//	filter[9*1+3*0+0],filter[9*1+3*0+1],filter[9*1+3*0+2],filter[9*2+3*0+0],filter[9*2+3*0+1],filter[9*2+3*0+2]);
			//printf("%d,%d,%d %d,%d,%d %d,%d,%d\n",filter[9*0+3*1+0],filter[9*0+3*1+1],filter[9*0+3*1+2],
			//	filter[9*1+3*1+0],filter[9*1+3*1+1],filter[9*1+3*1+2],filter[9*2+3*1+0],filter[9*2+3*1+1],filter[9*2+3*1+2]);
			//printf("%d,%d,%d %d,%d,%d %d,%d,%d\n\n",filter[9*0+3*2+0],filter[9*0+3*2+1],filter[9*0+3*2+2],
			//	filter[9*1+3*2+0],filter[9*1+3*2+1],filter[9*1+3*2+2],filter[9*2+3*2+0],filter[9*2+3*2+1],filter[9*2+3*2+2]);
			//ここから削除判定
			//周囲6近傍にボクセル無い個所があれば削除の可能性あり
			int hole_num = -1;	//はじめに見つけた穴のあいている箇所
			int x_num,y_num,z_num;
			for(int i=0;i<27;i++)
				if(filter[i]==-1){deleteFlag=true; hole_num=i; break;}
			if(deleteFlag){
				//注目ボクセルを除去しても穴があかないなら削除の可能性あり
				x_num=hole_num%3, y_num=(hole_num/3)%3, z_num=hole_num/9;
				erodeHoleCheck(filter,x_num,y_num,z_num);
				for(int i=0;i<27;i++)
					if(filter[i]==-1)return;
				//注目ボクセルを除去しても26近傍で連結しているなら削除
				int voxcel_num = -1;	//はじめに見つけたボクセルの存在する箇所
				for(int i=0;i<27;i++)
					if(filter[i]==2){voxcel_num=i; break;}
				x_num=voxcel_num%3, y_num=(voxcel_num/3)%3, z_num=voxcel_num/9;
				erodeConnectionCheck(filter,x_num,y_num,z_num);
				for(int i=0;i<27;i++)
					if(filter[i]==2)return;
				voxcel->state=0;
				voxcel->value=0;
			}
		}
	}else{
		if(voxcel->child==NULL)return;
		for(int i=0;i<voxcel->child_size;i++)
			erodeTVoxcel(voxcel->child[i],level,root_vox);
	}
}




/*----------------------------------------------------------------------------
	ボクセルを指定階層において膨張処理
	--------------------------------------------------------------------------
	引数levelで指定した階層数において、膨張処理を実行する
	この処理で新しく生成されたボクセルの状態は0となる（通常可視のボクセルは1）
----------------------------------------------------------------------------*/
//voxcel	ボクセル
//level		指定階層
//root_vox	ルートボクセル
void dilateTVoxcel(TVoxcel *voxcel, int level, TVoxcel *root_vox=NULL){
	if(root_vox==NULL)root_vox=voxcel;
	if(voxcel->level==level){
		if(voxcel->state>0){
			double vox_pos[3]; double vox_length[3];
			vox_pos[0]=(voxcel->max_x+voxcel->min_x)/2;
			vox_pos[1]=(voxcel->max_y+voxcel->min_y)/2;
			vox_pos[2]=(voxcel->max_z+voxcel->min_z)/2;
			vox_length[0]=voxcel->max_x-voxcel->min_x;
			vox_length[1]=voxcel->max_y-voxcel->min_y;
			vox_length[2]=voxcel->max_z-voxcel->min_z;
			for(int z=-1;z<=1;z++)
				for(int y=-1;y<=1;y++)
					for(int x=-1;x<=1;x++)
						devideTVoxcelByPoint(root_vox,
							vox_pos[0]+vox_length[0]*x,
							vox_pos[1]+vox_length[1]*y,
							vox_pos[2]+vox_length[2]*z,
							level,0);
		}
	}else{
		if(voxcel->child==NULL)return;
		for(int i=0;i<voxcel->child_size;i++)
			dilateTVoxcel(voxcel->child[i],level,root_vox);
	}
}

/*----------------------------------------------------------------------------
	ボクセルを指定階層において出力
	--------------------------------------------------------------------------
	引数levelで指定した階層数において、各ボクセルの中心座標を出力する
----------------------------------------------------------------------------*/
//voxcel	ボクセル
//level		指定階層
//fp		ファイルポインタ
void outputTVoxcel(TVoxcel *voxcel, int level, FILE *fp){
	if(voxcel->level==level){
		fprintf(fp,"%f %f %f\n",
			(voxcel->max_x+voxcel->min_x)/2,
			(voxcel->max_y+voxcel->min_y)/2,
			(voxcel->max_z+voxcel->min_z)/2);
	}else{
		if(voxcel->child==NULL)return;
		for(int i=0;i<voxcel->child_size;i++)
			outputTVoxcel(voxcel->child[i],level,fp);
	}
}

