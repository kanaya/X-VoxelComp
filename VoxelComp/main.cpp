// astyle --style=java -s2 -C -S -Y -p -H -k3 -W3 -y -j -c VoxelComp/main.cpp

/* 10年度卒業研究用プログラム『DTVoxcel(Divide Tree Voxcel)』の、CON(confidence)点群信頼度追加版、作成開始2011/02/07 */
/* DTVoxcel_REを流用(入力データはdouble型で自由)。OpenCVによるボクセルのXYZ軸断面の二次元画像を出す。 */
/* 対象オブジェクトの8分木は、XYZの内の最大の幅を利用し、空間分割する。REの信頼度は最深ボクセルのみ対応 */
/* 2011_02_10最終更新 */

#include <iostream>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <GLUT/glut.h>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <GL/glui.h>


//ボクセル系の関数
//ボクセル構造体の要素に、「内部に含まれる点の個数」「信頼度= (ボクセル内部の個数 / 全体の個数) * (８＾階層レベル)」

/*mask読み込み*/
// #include <windows.h>

#include <cmath>
#include <string>

//ボクセルの分割数(8分木構造体なら8で固定)
#define CHILD_VOXCEL_NUM 8


/*----------------------------------------------------------------------------
  ボクセル構造体
----------------------------------------------------------------------------*/
struct TVoxcel {
  public:
    int id;           //識別子
    TVoxcel *parent;      //親ノード（これがNullのとき、このノードは根）
    TVoxcel **child;      //子ノード（これがNullのとき、このノードは葉）
    int child_size;       //子ノードの数:8分木なら基本8
    int level;          //階層レベル（これが０のとき、このノードは根）
    double min_x, min_y, min_z; //この領域の最小の点
    double max_x, max_y, max_z; //この領域の最大の点
    int point_number;     //内包する点群の個数
    double value;       //このボクセルが持つ値(=階層数)
    double state;       //このボクセルの状態(-1:存在しない,0:処理中,1:存在,2:固定、など)
};

TVoxcel *root;            //八分木作成用ルートボクセル

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
    voxcel->parent = parent;
    voxcel->child = NULL;
    if (voxcel->child_size == 8) {  //8分木時：条件演算子→代入される物　＝　条件式　？　真：偽
      voxcel->min_x = (child_num % 2 == 0)     ? parent->min_x : (parent->min_x + parent->max_x) / 2;
      voxcel->min_y = (child_num / 2 % 2 == 0) ? parent->min_y : (parent->min_y + parent->max_y) / 2;
      voxcel->min_z = (child_num / 4 == 0)     ? parent->min_z : (parent->min_z + parent->max_z) / 2;
      voxcel->max_x = (child_num % 2 == 0)     ? (parent->min_x + parent->max_x) / 2 : parent->max_x;
      voxcel->max_y = (child_num / 2 % 2 == 0) ? (parent->min_y + parent->max_y) / 2 : parent->max_y;
      voxcel->max_z = (child_num / 4 == 0)     ? (parent->min_z + parent->max_z) / 2 : parent->max_z;
    }
    else { //八分木ではない時
      voxcel->min_x = 0, voxcel->min_y = 0, voxcel->min_z = 0;
      voxcel->max_x = 1, voxcel->max_y = 1, voxcel->max_z = 1;
    }
    voxcel->point_number = 0;
    voxcel->value = 0;
    voxcel->state = state;
  }
  else {  //ルートボクセル用
    voxcel->level = 0;
    voxcel->id    = 0;
    voxcel->parent = NULL;
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
  for (int i = 0; i < voxcel->child_size; i++) { //１〜８(8分木なら8、4分木なら4が入る)
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
  if (voxcel->child != 0) {
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
  TVoxcel *tmp_vox = voxcel;  // 現在返す候補のボクセル
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
  if (includePointInTVoxcel(voxcel, x, y, z)) {
    if (level > voxcel->level) {
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
//voxcel    ボクセル
//pre_state   置き換えるボクセルの状態値
//set_state   変更後の状態の値
void replacedTVoxcelState(TVoxcel *voxcel, double pre_state, double set_state) {
  if (voxcel->state == pre_state) {
    voxcel->state = set_state;
  }
  if (voxcel->child != NULL)
    for (int i = 0; i < voxcel->child_size; i++) {
      replacedTVoxcelState(voxcel->child[i], pre_state, set_state);
    }
}

/*----------------------------------------------------------------------------
  ボクセルを削除(サブルーチン)
  --------------------------------------------------------------------------
  ボクセルを子要素を消し，自分も消す．
----------------------------------------------------------------------------*/
//voxcel: 削除するボクセル：cppからはdeleteを経由してrootが送られる
void Rec_deleteTVoxcel(TVoxcel *voxcel) {
  if (voxcel->child != NULL) {
    for (int i = 0; i < 8; i++) {
      //TVoxcel* voxcel = voxcel->child[i];
      //Rec_deleteTVoxcel(voxcel);  // 子以降を削除します
      if (voxcel->child[i] != NULL) {
        Rec_deleteTVoxcel(voxcel->child[i]);  // 子以降を削除します
      }
    }
  }
  if (voxcel->id < 80) {
    fprintf(stderr, "->");  //確認用
  }
  delete voxcel;  // 最後に自分を削除
}


/*----------------------------------------------------------------------------
  縮退処理において注目ボクセルを削除すると26近傍で穴があくかのチェック
  --------------------------------------------------------------------------
  引数filterの配列に対し、穴のあいた箇所を指定の条件ですべて通ることができるか
  通ることができるのは値が2の要素で、なおかつ注目ボクセルや端の8箇所のボクセル
  を除く箇所である。すべて通るならば、値が-1の箇所はすべて0に書き換えられる
----------------------------------------------------------------------------*/
//filter  フィルタ
//x,y,z   現在の位置
void erodeHoleCheck(int *filter, int x, int y, int z) {
  if (filter[z * 9 + y * 3 + x] == -1) {
    filter[z * 9 + y * 3 + x] = 0;
  }
  else if (filter[z * 9 + y * 3 + x] == 2) {
    filter[z * 9 + y * 3 + x] = 3;
  }
  else {
    return;
  }
  for (int k = -1; k <= 1; k++)
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
        if (((z + k) * 9 + (y + j) * 3 + (x + i)) >= 0 && ((z + k) * 9 + (y + j) * 3 + (x + i)) < 27 && !(
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 0 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 2 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 6 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 8 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 18 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 20 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 24 ||
              ((z + k) * 9 + (y + j) * 3 + (x + i)) == 26) ) {
          erodeHoleCheck(filter, x + i, y + j, z + k);
        }
}

/*----------------------------------------------------------------------------
  縮退処理において注目ボクセルを削除しても26近傍で連結しているかのチェック
  --------------------------------------------------------------------------
  引数filterの配列に対し、値が2か3の要素を26近傍の移動ですべて通ることができるか
  すべて通るならば、値が2か3の箇所はすべて4に書き換えられる
----------------------------------------------------------------------------*/
//filter  フィルタ
//x,y,z   現在の位置
void erodeConnectionCheck(int *filter, int x, int y, int z) {
  if (filter[z * 9 + y * 3 + x] == 2 || filter[z * 9 + y * 3 + x] == 3) {
    filter[z * 9 + y * 3 + x] = 4;
  }
  else {
    return;
  }
  for (int k = -1; k <= 1; k++)
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
        if (((z + k) * 9 + (y + j) * 3 + (x + i)) >= 0 && ((z + k) * 9 + (y + j) * 3 + (x + i)) < 27) {
          erodeConnectionCheck(filter, x + i, y + j, z + k);
        }
}

/*----------------------------------------------------------------------------
  ボクセルを指定階層において縮退処理
  --------------------------------------------------------------------------
  引数levelで指定した階層数において、縮退処理を実行する
  この処理で新しく生成されたボクセルの状態は0となる（通常可視のボクセルは1）
----------------------------------------------------------------------------*/
//voxcel  ボクセル
//level   指定階層
//root_vox  ルートボクセル
void erodeTVoxcel(TVoxcel *voxcel, int level, TVoxcel *root_vox = NULL) {
  if (root_vox == NULL) {
    root_vox = voxcel;
  }
  if (voxcel->level == level) {
    if (voxcel->state == 1) {
      int filter[27];
      TVoxcel *tmp_vox;
      bool deleteFlag = false;
      double vox_pos[3];
      double vox_length[3];
      vox_pos[0] = (voxcel->max_x + voxcel->min_x) / 2;
      vox_pos[1] = (voxcel->max_y + voxcel->min_y) / 2;
      vox_pos[2] = (voxcel->max_z + voxcel->min_z) / 2;
      vox_length[0] = voxcel->max_x - voxcel->min_x;
      vox_length[1] = voxcel->max_y - voxcel->min_y;
      vox_length[2] = voxcel->max_z - voxcel->min_z;
      //printf("%f,%f,%f(%f,%f,%f)\n",vox_pos[0],vox_pos[1],vox_pos[2],vox_length[0],vox_length[1],vox_length[2]);
      //注目ボクセルの周囲のボクセル状況を記録
      //-1:6近傍のボクセルなし,0:26近傍のボクセルなしか注目ボクセル
      //1:ボクセル不可視,2:ボクセル存在
      for (int z = -1; z <= 1; z++) {
        for (int y = -1; y <= 1; y++) {
          for (int x = -1; x <= 1; x++) {
            tmp_vox = getTVoxcel(root_vox,
                                 vox_pos[0] + vox_length[0] * x,
                                 vox_pos[1] + vox_length[1] * y,
                                 vox_pos[2] + vox_length[2] * z,
                                 level);
            if (tmp_vox != 0
                && (9 * (z + 1) + 3 * (y + 1) + (x + 1)) != 13
                && tmp_vox->state != -1) {
              if (tmp_vox->state == 0) {
                filter[9 * (z + 1) + 3 * (y + 1) + (x + 1)] = 1;
              }
              else {
                filter[9 * (z + 1) + 3 * (y + 1) + (x + 1)] = 2;
              }
            }
            else {
              //printf("%d\n",(9*(z+1)+3*(y+1)+(x+1)));
              if ( (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 4 ||
                   (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 10 ||
                   (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 12 ||
                   (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 14 ||
                   (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 16 ||
                   (9 * (z + 1) + 3 * (y + 1) + (x + 1)) == 22) {
                filter[9 * (z + 1) + 3 * (y + 1) + (x + 1)] = -1;
              }
              else {
                filter[9 * (z + 1) + 3 * (y + 1) + (x + 1)] = 0;
              }
            }
          }
        }
      }
      //printf("%d,%d,%d %d,%d,%d %d,%d,%d\n",filter[9*0+3*0+0],filter[9*0+3*0+1],filter[9*0+3*0+2],
      //  filter[9*1+3*0+0],filter[9*1+3*0+1],filter[9*1+3*0+2],filter[9*2+3*0+0],filter[9*2+3*0+1],filter[9*2+3*0+2]);
      //printf("%d,%d,%d %d,%d,%d %d,%d,%d\n",filter[9*0+3*1+0],filter[9*0+3*1+1],filter[9*0+3*1+2],
      //  filter[9*1+3*1+0],filter[9*1+3*1+1],filter[9*1+3*1+2],filter[9*2+3*1+0],filter[9*2+3*1+1],filter[9*2+3*1+2]);
      //printf("%d,%d,%d %d,%d,%d %d,%d,%d\n\n",filter[9*0+3*2+0],filter[9*0+3*2+1],filter[9*0+3*2+2],
      //  filter[9*1+3*2+0],filter[9*1+3*2+1],filter[9*1+3*2+2],filter[9*2+3*2+0],filter[9*2+3*2+1],filter[9*2+3*2+2]);
      //ここから削除判定
      //周囲6近傍にボクセル無い個所があれば削除の可能性あり
      int hole_num = -1;  //はじめに見つけた穴のあいている箇所
      int x_num, y_num, z_num;
      for (int i = 0; i < 27; i++) {
        if (filter[i] == -1) {
          deleteFlag = true;
          hole_num = i;
          break;
        }
      }
      if (deleteFlag) {
        //注目ボクセルを除去しても穴があかないなら削除の可能性あり
        x_num = hole_num % 3, y_num = (hole_num / 3) % 3, z_num = hole_num / 9;
        erodeHoleCheck(filter, x_num, y_num, z_num);
        for (int i = 0; i < 27; i++) {
          if (filter[i] == -1) {
            return;
          }
        }
        //注目ボクセルを除去しても26近傍で連結しているなら削除
        int voxcel_num = -1;  //はじめに見つけたボクセルの存在する箇所
        for (int i = 0; i < 27; i++) {
          if (filter[i] == 2) {
            voxcel_num = i;
            break;
          }
        }
        x_num = voxcel_num % 3;
        y_num = (voxcel_num / 3) % 3;
        z_num = voxcel_num / 9;
        erodeConnectionCheck(filter, x_num, y_num, z_num);
        for (int i = 0; i < 27; i++)
          if (filter[i] == 2) {
            return;
          }
        voxcel->state = 0;
        voxcel->value = 0;
      }
    }
  }
  else {
    if (voxcel->child == NULL) {
      return;
    }
    for (int i = 0; i < voxcel->child_size; i++) {
      erodeTVoxcel(voxcel->child[i], level, root_vox);
    }
  }
}




/*----------------------------------------------------------------------------
  ボクセルを指定階層において膨張処理
  --------------------------------------------------------------------------
  引数levelで指定した階層数において、膨張処理を実行する
  この処理で新しく生成されたボクセルの状態は0となる（通常可視のボクセルは1）
----------------------------------------------------------------------------*/
//voxcel  ボクセル
//level   指定階層
//root_vox  ルートボクセル
void dilateTVoxcel(TVoxcel *voxcel, int level, TVoxcel *root_vox = NULL) {
  if (root_vox == NULL) {
    root_vox = voxcel;
  }
  if (voxcel->level == level) {
    if (voxcel->state > 0) {
      double vox_pos[3];
      double vox_length[3];
      vox_pos[0] = (voxcel->max_x + voxcel->min_x) / 2;
      vox_pos[1] = (voxcel->max_y + voxcel->min_y) / 2;
      vox_pos[2] = (voxcel->max_z + voxcel->min_z) / 2;
      vox_length[0] = voxcel->max_x - voxcel->min_x;
      vox_length[1] = voxcel->max_y - voxcel->min_y;
      vox_length[2] = voxcel->max_z - voxcel->min_z;
      for (int z = -1; z <= 1; z++)
        for (int y = -1; y <= 1; y++)
          for (int x = -1; x <= 1; x++)
            devideTVoxcelByPoint(root_vox,
                                 vox_pos[0] + vox_length[0]*x,
                                 vox_pos[1] + vox_length[1]*y,
                                 vox_pos[2] + vox_length[2]*z,
                                 level, 0);
    }
  }
  else {
    if (voxcel->child == NULL) {
      return;
    }
    for (int i = 0; i < voxcel->child_size; i++) {
      dilateTVoxcel(voxcel->child[i], level, root_vox);
    }
  }
}

/*----------------------------------------------------------------------------
  ボクセルを指定階層において出力
  --------------------------------------------------------------------------
  引数levelで指定した階層数において、各ボクセルの中心座標を出力する
----------------------------------------------------------------------------*/
//voxcel  ボクセル
//level   指定階層
//fp    ファイルポインタ
void outputTVoxcel(TVoxcel *voxcel, int level, FILE *fp) {
  if (voxcel->level == level) {
    fprintf(fp, "%f %f %f\n",
            (voxcel->max_x + voxcel->min_x) / 2,
            (voxcel->max_y + voxcel->min_y) / 2,
            (voxcel->max_z + voxcel->min_z) / 2);
  }
  else {
    if (voxcel->child == NULL) {
      return;
    }
    for (int i = 0; i < voxcel->child_size; i++) {
      outputTVoxcel(voxcel->child[i], level, fp);
    }
  }
}




// #define DO_WHEN_DEBUG(x) (x)
#define DO_WHEN_DEBUG(x) ((void)0)

#define MAX_DEVIDE_FREQUENCY 9 //八分木最大(MAX)分割(DIVIDE)回数(FREQUENCY)
#define SCAN3D_POINT_DATA "/Users/kanaya/Documents/VoxelComp/VoxelComp/f_church.txt"
#define OUTPUTDATA "./f_church_processed.txt"


const int nil = 0;

/*----------------------
 変数定義
 ----------------------*/
/****[GL定義]****/
const char *gl_title = "Voxcel_View_GL";    //GLウィンドウのタイトル
int gl_width   = 640;         //GLウィンドウの幅
int gl_height  = 800;         //GLウィンドウの高さ
double vox_min[3], vox_max[3];      //対象オブジェクトの領域の最小(x,y,z)と最大(x,y,z)
double vox_max_divide;          //対象オブジェクトの領域の最大の幅、正規格子にする為に使用
static std::vector<double> data_3d;   //x,y,zの順に点データ格納
GLuint dispList;            //ディスプレイリスト
double dist = 100;            //視点までの距離:表示される前に一応変更される。
int vox_value = 5;            //描画するボクセルの階層(8^-vox_valueの大きさ) [初期設定]：キーa/Aで変更
double rotX = 0, rotY = 0, rotZ = 0;      //回転量[初期設定]：各キーxyz/XYZで変更、360度
double XX, YY, ZZ;            //X・Y・Z軸の断面における『中心』。
double DX;                //断面幅：断面を取る幅対象オブジェクトの最大辺を一辺とする空間(の１辺)を分割したもの．
double DP, DM;              //断面幅の高い値/低い値（暫定）:フラッグによりX/Y/Zのそれぞれの幅に変更される
double glx, gly, glz;           //CVから点を追加する際に、一時的に指定した値をいれておく変数。
//GL:count
int countDX, countDY, countDZ;      //任意断面の移動数管理：a/Aの整合性を取る。細かくした場合、値が小さい方の幅に移る。
int confidence_p = 0;           //点群の粗密から信頼度を計算する時にしようする。
int confidence_max = 0;         //点群の粗密の内、選択中の階層ボクセルの中で最大の個数
//GL:Flag
int dVFlag = 0;             //・g/G:drawVoxcelのフラッグ。0なら信頼度無し、1なら信頼度込み。
int acFlag = 0;             //・l/L:Voxcel全体を表示するのか、任意断面のみ表示するか。all ot choise Flag:未実装
int fcount = 0;             //グローバルでとにかく確認する為のもの。


/****[CV定義]****/
const char *cv_title = "Voxcel_View_CV";
int cv_width = 1024;                  //CVウィンドウの幅
int cv_height = 1024;                 //CVウィンドウの高さ
CvSize window_size = {cv_width, cv_height};       //窓サイズ
static int CVcnt = 0;                     //CV反映・CvPoint型データ・カウント用
int CV3dcnt = 0;                      //CV反映・double型データ・カウント用。
const int dimension = MAX_DEVIDE_FREQUENCY;       //分割数、初期のデフォルトは10．
static int DD = 0;                    //採取する断面の種類読み分け用[初期設定]：0:X、1:Y、2:Z,data_3d[DD]
static int DDX = 1, DDY = 2;                  //同上[初期設定]：X,Y,Z断面時。X軸採取ならYとZなので3次元data_3d[i+DDX or i+DDY]。
IplImage *imgA;                     //CV用IPLイメージを作成
CvPoint *CVPts;                     //OpenCV用の二次元座標群、CVPts[cnt].xとの形、int型の値。
//CV:Flag
int DeleteFlag = 0;                   //・m/M:CVのﾏｳｽｸﾘｯｸ時、0:点追加.1:点削除
int PlusFlag = 0;                   //・n/N:CVで点追加時、0：指定断面の中央、1：指定ボクセル内の全最深ボクセルの中央
int PixelFlag = 0;                    //・o/O:CVでピクセル：0：表示する。１：表示しない。
int ConFlag = 0;                      //・p/P:CVで信頼度、0：表示する。１：表示しない。
int ReviewFlag = 0;                   //・q/Q：CVでクリックの旅に再表示するか：0:表示する。1:しない。

/******[GLUI定義]******/
GLUI *control;
float rotate[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
};
int window_id;
float gluiXYZ[3] = { 0.0, 0.0, 0.0 };


/***[プロトタイプ宣言]***/
void drawVoxcel0(TVoxcel *voxcel, int value);           //GLのボクセル描画・カラー：基本RGB
void drawVoxcel1(TVoxcel *voxcel, int value);           //GLのボクセル描画・カラー：点群粗密信頼度HSV
void value_counter_con(void);                   //GLの点群信頼度算出用
void display(void);                         //GLの表示
void counterCpoint(TVoxcel *voxcel, int value);           //GLの点群信頼度、MAX算出用
void reshape(int x, int y);                     //画面再描画
void timer();                           //timer関数：常時Reディスプレイ
void DRAW_TREE_4( IplImage *clone, CvPoint *KData);//CV:四分木：統括
void DRAW_TREE_4c2( IplImage *imgA, CvPoint *pts , int x, int y, int dim); //CV:四分木信頼度
void DRAW_TREE_4c1( IplImage *imgA, CvPoint *pts , int x, int y, int dim); //CV:四分木信頼度
void DRAW_TREE_4v( IplImage *imgA, CvPoint *pts , int x, int y, int dim); //CV:四分木ボクセル
void DRAW_TREE_4p( IplImage *imgA, CvPoint *pts , int x, int y, int dim); //CV:四分木ポイント
void DRAW_TREE_4l( IplImage *imgA, CvPoint *pts , int x, int y, int dim); //CV:四分木ライン
void KnownPoint( void *imgA, CvPoint *KData );            //CV：既知点表示
CvPoint *CVset(void);                       //CV用の座標点、3次元点群からセット
static void MOUSE(int event, int x, int y, int flags, void *imgA);  //マウス
static void keyboard(unsigned char key, int x, int y);        //キーボード
static void DAxis(int Flag);                    //断面デフォルト
void gluiCallbackExit(int num);                   //GLUI:EXIT
void gluiCallbackDef(int num);                    //GLUI:再描画
void FileOutput(void);                        //ファイル出力
void ReFileInput(const char *Input_data);                 //ファイル入力
void DataRead(void);                        //データ閲覧
void Make_Voxcel(void);                       //ボクセルデータ作成
void point_plus(int plus_x, int plus_y, void *imgA);        //CVで点を追加した際に実行される関数
void read_Point(const char *Input_data);                  //点群ファイルから、データ型に移行[初期設定]




/*---------------------------------------------------
 メインメソッド
 -----------------------------------------------------*/
int main(int argc, char **argv) {
  fprintf(stderr, "GLWindowSize=(%d, %d)\n", gl_width, gl_height);
  fprintf(stderr, "CVWindowSize=(%d, %d)\n", cv_width, cv_height);

  /*** OpenGL窓作成 ***/
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(gl_width, gl_height);
  glutInitWindowPosition(100, 100);

  //GLUI改変スタート
  //  glutCreateWindow(gl_title);
  window_id = glutCreateWindow( gl_title );//GLUI試験中
  GLUI_Master.set_glutIdleFunc( nil );//GLUI試験中

  /*** OpenCV窓作成 ***/
  imgA = cvCreateImage(window_size, IPL_DEPTH_8U, 3);
  cvSet(imgA, cvScalarAll (255), 0);
  cvNamedWindow(cv_title, CV_WINDOW_AUTOSIZE);

  /*** ファイルから点の入力 ***/
  read_Point(SCAN3D_POINT_DATA);

  /*** 点データからボクセルの作成 ***/
  Make_Voxcel();
  dist = (vox_max[2] - vox_min[2]) * 2; //カメラ位置変更


  /*** CV用断面の値設定 ***///正規格子変更考慮箇所
  DX = vox_max_divide / (pow(2, (double)vox_value));  //断面幅[初期設定]
  XX = vox_min[0] + DX / 2;                   //X軸断面[初期設定]：最も小さいX値のボクセル中心に移動
  YY = vox_min[1] + DX / 2;                   //Y軸断面[初期設定]：最も小さいY値のボクセル中心に移動
  ZZ = vox_min[2] + DX / 2;                   //Z軸断面[初期設定]：最も小さいZ値のボクセル中心に移動
  DM = XX - DX / 2;                       //断面の小さい値[初期設定]：
  DP = XX + DX / 2;                       //断面の大きい値[初期設定]：
  countDX = countDY = countDZ = 0;                  //断面幅の移動数
  fprintf(stderr, "\n初期断面幅[%f],初期X軸断面中央[%f],初期Y軸断面中央[%f],初期Z軸断面中央[%f]\n", DX, XX, YY, ZZ);
  fprintf(stderr, "表示階層数[%d]:断面幅[%f]\n", vox_value, DX);

  //デプステスト:多角形の前後関係を把握する
  glClearColor(1, 1, 1, 1.0);
  glClear(GL_DEPTH_BUFFER_BIT);                     //デプスバッファのクリア
  glEnable(GL_DEPTH_TEST);                        //デプステストを可能に
  glEnable(GL_NORMALIZE);                         //法線を自動的に正規化

  //スムーズシェーディング
  glShadeModel(GL_SMOOTH);

  //関数の初期設定
  GLUI_Master.set_glutIdleFunc(timer);
  glutDisplayFunc(display);
  control = GLUI_Master.create_glui_subwindow(window_id, GLUI_SUBWINDOW_BOTTOM);//GLUI実験中
  GLUI_Panel *transPanel = control->add_panel( "Translation Control");
  control->add_column_to_panel( transPanel, false );
  control->add_rotation_to_panel(transPanel, "Rotation", rotate); //GLUI:rotation
  control->add_button("Exit", 0, gluiCallbackExit);//GLUI:
  control->add_button("Default", 0, gluiCallbackDef);//GLUI:name,id,
  control->add_column_to_panel( transPanel, false );
  control->add_translation_to_panel( transPanel, "Translation X-Y", GLUI_TRANSLATION_XY, gluiXYZ );
  control->add_column_to_panel( transPanel, false );
  control->add_translation_to_panel( transPanel, "Translation Z", GLUI_TRANSLATION_Z, &gluiXYZ[2]);
  control->set_main_gfx_window( window_id);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  cvSetMouseCallback(cv_title, MOUSE, imgA);
  MOUSE(CV_EVENT_LBUTTONUP, 0, 0, 0, imgA);  // TEST TEST TEST
  glutPostRedisplay();
  glutMainLoop();

  cvReleaseImage(& imgA);
  cvDestroyWindow(cv_title);
  return 0;
}

/*----------------------
 ボクセル描画：基本
 ----------------------*/
//voxcel: トップレベルボクセル
//value:  描画するボクセルの値
void drawVoxcel0(TVoxcel *voxcel, int value) { // 再帰
  int Con = 0;  // ワイヤー[0]orソリッド[1]
  int choise = 0;  // 選択されるXYZ断面かどうか：１なら断面候補
  double color_x, color_y, color_z;
  color_x = color_y = color_z = 0;
  if (voxcel->value == value) {
    glPushMatrix();
    glColor4f(1.0, 1.0, 1.0, 0.0);
    Con = 0;
    if (voxcel->min_x <= XX && XX < voxcel->max_x) {
      glColor4f(color_x += 0.9 * (voxcel->max_z - vox_min[2]) / (vox_max[2] - vox_min[2]) + 0.1,
                color_y,
                color_z,
                0.5);
      choise = 1;//カラー設定X断面
    }
    if (voxcel->min_y <= YY && YY < voxcel->max_y) {
      glColor4f(color_x,
                color_y += 0.9 * (voxcel->max_z - vox_min[2]) / (vox_max[2] - vox_min[2]) + 0.1,
                color_z,
                0.5);
      choise = 1;//カラー設定Y断面
    }
    if (voxcel->min_z <= ZZ && ZZ < voxcel->max_z) {
      glColor4f(color_x,
                color_y,
                color_z += 0.9 * (voxcel->max_z - vox_min[2]) / (vox_max[2] - vox_min[2]) + 0.1,
                0.8);
      choise = 1;//カラー設定Z断面
    }
    // ここまでで色と、ワイヤーorソリッドを決める。
    if (DD == 0) {
      if (voxcel->min_x <= XX && XX < voxcel->max_x) {
        Con = 1;
      }
    }
    if (DD == 1) {
      if (voxcel->min_y <= YY && YY < voxcel->max_y) {
        Con = 1;
      }
    }
    if (DD == 2) {
      if (voxcel->min_z <= ZZ && ZZ < voxcel->max_z) {
        Con = 1;
      }
    }
    glTranslatef((voxcel->max_x + voxcel->min_x) / 2,
                 (voxcel->max_y + voxcel->min_y) / 2,
                 (voxcel->max_z + voxcel->min_z) / 2);  // 移動して
    glScalef(voxcel->max_x - voxcel->min_x,
             voxcel->max_y - voxcel->min_y,
             voxcel->max_z - voxcel->min_z);  // 大きさ指定して
    if (acFlag == 0) {
      if (Con == 0) {
        glutWireCube(1);        //ワイヤーを描く
      }
    }
    if (acFlag == 1) {
      if (Con == 0) {
        if (choise == 1) {
          glutWireCube(1);  // 選択断面はその他のワイヤーが消えても描く
        }
      }
    }
    if (Con == 1) {
      glutSolidCube(1);             //ソリッドを描く
    }
    glPopMatrix();
  }
  else {   // 再帰：ボクセルのvalueがvalue（指定した階層）ではなく、さらに〜分木の構成の場合、
    if (voxcel->child != nil) { //N分木ならN回、自分の子を行う。例えば5階層のボクセル一覧を表示するなら、5階層までもぐって表示させる。
      for (int i = 0; i < voxcel->child_size; i++) {
        drawVoxcel0(voxcel->child[i], value);
      }
    }
  }
}


/*----------------------
 ボクセル描画改良1：点群による信頼度
 ----------------------*/
void drawVoxcel1(TVoxcel *voxcel, int value) {
  int Con = 0; //ワイヤーorソリッド
  int choise = 0; //選択されるXYZ断面かどうか：１なら断面候補
  int HSV_Hi;               //HSV色空間変換用
  double HSV_H, HSV_f, HSV_p, HSV_q, HSV_t; //HSV色空間変換用
  double siken = 0;
  if (voxcel->value == value) {
    confidence_p = 0;
    //siken=((double)voxcel->point_number/(double)confidence_max+0.8)*((double)voxcel->point_number/(double)confidence_max+0.8)/3.24;
    siken = (double)voxcel->point_number / (double)confidence_max;
    HSV_H = siken * 270; //ボクセル内包個数/指定階層ボクセル内包最大個数を、０°〜270°で表現する
    HSV_Hi = (int)(HSV_H / 60) % 6; //HSV色の360度中270度で表現する。
    HSV_f = ((double)HSV_H / 60) - (double)HSV_Hi;
    HSV_p = 0;
    HSV_q = (1 - (double)HSV_f);
    HSV_t = (1 - (1 - (double)HSV_f));
    glPushMatrix();
    if ((voxcel->min_x <= XX && XX < voxcel->max_x)
        || (voxcel->min_y <= YY && YY < voxcel->max_y)
        || (voxcel->min_z <= ZZ && ZZ < voxcel->max_z)) { //カラー設定断面:HSV表記
      if (HSV_Hi == 0) {
        glColor4f(1.0, HSV_t, HSV_p, 1.0);
      }
      if (HSV_Hi == 1) {
        glColor4f(HSV_q, 1.0, HSV_p, 1.0);
      }
      if (HSV_Hi == 2) {
        glColor4f(HSV_p, 1.0, HSV_t, 1.0);
      }
      if (HSV_Hi == 3) {
        glColor4f(HSV_p, HSV_q, 1.0, 1.0);
      }
      if (HSV_Hi == 4) {
        glColor4f(HSV_t, HSV_p, 1.0, 1.0);
      }
      if (HSV_Hi == 5) {
        glColor4f(1.0, HSV_p, HSV_q, 1.0);
      }
      choise = 1; //断面候補
    }
    else {
      glColor4f(1.0, 1.0, 1.0, 0.0), Con = 0; //ここまでで色と、ワイヤーorソリッドを決める。
    }
    if (DD == 0)if (voxcel->min_x <= XX && XX < voxcel->max_x) {
        Con = 1;
      }
    if (DD == 1)if (voxcel->min_y <= YY && YY < voxcel->max_y) {
        Con = 1;
      }
    if (DD == 2)if (voxcel->min_z <= ZZ && ZZ < voxcel->max_z) {
        Con = 1;
      }
    glTranslatef((voxcel->max_x + voxcel->min_x) / 2, (voxcel->max_y + voxcel->min_y) / 2, (voxcel->max_z + voxcel->min_z) / 2); //移動して
    glScalef(voxcel->max_x - voxcel->min_x, voxcel->max_y - voxcel->min_y, voxcel->max_z - voxcel->min_z); //大きさ指定して
    if (acFlag == 0)if (Con == 0) {
        glutWireCube(1);  //ワイヤーを描く
      }
    if (acFlag == 1)if (Con == 0)if (choise == 1) {
          glutWireCube(1);  //選択断面はその他のワイヤーが消えても描く
        }
    if (Con == 1) {
      glutSolidCube(1);  //ソリッドを描く
    }
    glPopMatrix();
  }
  else {   //再帰：ボクセルのvalueがvalue（指定した階層）ではなく、さらに〜分木の構成の場合、
    if (voxcel->child != nil) { //N分木ならN回、自分の子を行う。例えば5階層のボクセル一覧を表示するなら、5階層までもぐって表示させる。
      for (int i = 0; i < voxcel->child_size; i++) {
        drawVoxcel1(voxcel->child[i], value);
      }
    }
  }
}



/*----------------------
 画面描画
 ----------------------*/
void display(void) {
  //色の設定
  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //オブジェクトの描画
  glDisable(GL_LIGHTING);
  glViewport(0, 0, gl_width, gl_height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40, gl_width / (GLdouble)gl_height, 0.01, 300000);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  //カメラの視点設定
  gluLookAt((vox_min[0] + vox_max[0]) / 2 - gluiXYZ[0], (vox_min[1] + vox_max[1]) / 2 - gluiXYZ[1], vox_max[2] + dist - gluiXYZ[2],
            (vox_min[0] + vox_max[0]) / 2 - gluiXYZ[0], (vox_min[1] + vox_max[1]) / 2 - gluiXYZ[1], (vox_min[2] + vox_max[2]) / 2, 0, 1, 0); //位置、注視中央、上方向

  //カメラ角度変更
  glTranslatef((vox_min[0] + vox_max[0]) / 2, (vox_min[1] + vox_max[1]) / 2, (vox_min[2] + vox_max[2]) / 2);
  glRotatef(rotX, -1, 0, 0);
  glRotatef(rotY, 0, -1, 0);
  glRotatef(rotZ, 0, 0, -1);
  glTranslatef(-(vox_min[0] + vox_max[0]) / 2, -(vox_min[1] + vox_max[1]) / 2, -(vox_min[2] + vox_max[2]) / 2);

  //GLUI
  glPushMatrix();
  glTranslatef((vox_min[0] + vox_max[0]) / 2, (vox_min[1] + vox_max[1]) / 2, (vox_min[2] + vox_max[2]) / 2);
  glMultMatrixf(rotate);//GLUI
  glTranslatef(-(vox_min[0] + vox_max[0]) / 2, -(vox_min[1] + vox_max[1]) / 2, -(vox_min[2] + vox_max[2]) / 2);

  //ボクセル描画(ヘッダー)
  if (dVFlag == 0) {
    drawVoxcel0(root, vox_value);
  }
  if (dVFlag == 1) {
    drawVoxcel1(root, vox_value);
  }

  glPopMatrix();


  glutSwapBuffers();
}


/*----------------------
 画面再描画
 ----------------------*/
void reshape(int x, int y) {
  gl_width = x, gl_height = y;
  glutPostRedisplay();
}

/*----------------------
 timer関数
 ----------------------*/
void timer() {
  glutPostRedisplay();
}

/*-------------------------------------------------------------
 [OpenCV用関数]
 -------------------------------------------------------------*/

/*****************[CV表示：4分木:統括]**************************/
void DRAW_TREE_4(IplImage *clone, CvPoint *KData) {
  // 点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。
  CvPoint *pts = (CvPoint *)cvAlloc(sizeof(CvPoint) * 200);
  int i = 0;

  fprintf(stderr, "4分木処理実行開始\n"); // 確認用テンプレ
  DO_WHEN_DEBUG(fprintf(stderr, "1\n"));
  // 信頼度低
  if (ConFlag == 0) {
    DO_WHEN_DEBUG(fprintf(stderr, "CVcnt == %d >>>", CVcnt));
    for (i = 0; i < CVcnt; i++) { // CVcntはKdataの行数。つまり読み取った点群の行数
      DO_WHEN_DEBUG(fprintf(stderr, "."));
      pts[0] = cvPoint(0, 0);
      pts[1] = cvPoint(clone->width, clone->height);
      if (KData[i].x != nil) {
        DRAW_TREE_4c2(clone, pts, KData[i].x, KData[i].y, vox_value);
      }
      if (i % 100 == 0) {
        fprintf(stderr, "\r4分木：信頼度『低』実行[%d]行目：x[%d]y[%d] ", i, KData[i].x, KData[i].y);  // 確認用テンプレ
      }
    }
    // fprintf(stderr, ":完了[%d]行：x[%d]y[%d] ", i, KData[i - 1].x, KData[i - 1].y); // 確認用テンプレ
  }
  else {
    fprintf(stderr, "\n");
  }
  //信頼度高
  DO_WHEN_DEBUG(fprintf(stderr, "2\n"));
  if (ConFlag == 0) {
    for (i = 0; i < CVcnt; i++) { //CVcntはKdataの行数。つまり読み取った点群の行数
      pts[0] = cvPoint(0, 0);
      pts[1] = cvPoint(clone->width, clone->height);
      if (KData[i].x != nil) {
        DRAW_TREE_4c1( clone, pts, KData[i].x, KData[i].y, vox_value);
      }
      if (i % 100 == 0) {
        fprintf(stderr, "\r4分木：信頼度『高』実行[%d]行目：x[%d]y[%d] ", i, KData[i].x, KData[i].y);  //確認用テンプレ
      }
    }
    // fprintf(stderr, ":完了[%d]行：x[%d]y[%d] ",i,KData[i-1].x,KData[i-1].y);//確認用テンプレ
  }
  else {
    fprintf(stderr, "\n");
  }
  //ピクセル
  DO_WHEN_DEBUG(fprintf(stderr, "3\n"));
  if (PixelFlag == 0) { //見難い ON/OFF
    for (i = 0; i < CVcnt; i++) { // CVcntはKdataの行数。つまり読み取った点群の行数
      pts[0] = cvPoint(0, 0);
      pts[1] = cvPoint(clone->width, clone->height);
      if (KData[i].x != nil) {
        DRAW_TREE_4v( clone, pts, KData[i].x, KData[i].y, vox_value);
      }
      if (i % 100 == 0) {
        fprintf(stderr, "\r4分木：ボクセル着色実行[%d]行目：x[%d]y[%d] ", i, KData[i].x, KData[i].y);  //確認用テンプレ
      }
    }
    // fprintf(stderr, ":完了[%d]行：x[%d]y[%d] ",i,KData[i-1].x,KData[i-1].y);//確認用テンプレ
  }
  else {
    fprintf(stderr, "\n");
  }
  //点群
  DO_WHEN_DEBUG(fprintf(stderr, "4\n"));
  for (i = 0; i < CVcnt; i++) { //CVcntはKdataの行数。つまり読み取った点群の行数
    pts[0] = cvPoint(0, 0);
    pts[1] = cvPoint(clone->width, clone->height);
    if (KData[i].x != nil) {
      DRAW_TREE_4p( clone, pts, KData[i].x, KData[i].y, MAX_DEVIDE_FREQUENCY);
    }
    if (i % 100 == 0) {
      fprintf(stderr, "\r4分木：点群実行[%d]行目：x[%d]y[%d] ", i, KData[i].x, KData[i].y);  //確認用テンプレ
    }
  }
  // fprintf(stderr, ":完了[%d]行：x[%d]y[%d]\n",i,KData[i-1].x,KData[i-1].y);
  //ライン
  DO_WHEN_DEBUG(fprintf(stderr, "5\n"));
  for ( i = 0; i < CVcnt; i++) { //CVcntはKdataの行数。つまり読み取った点群の行数
    pts[0] = cvPoint(0, 0);
    pts[1] = cvPoint(clone->width, clone->height);
    if (KData[i].x != nil) {
      DRAW_TREE_4l( clone, pts, KData[i].x, KData[i].y, vox_value);
    }
    if (i % 100 == 0) {
      fprintf(stderr, "\r4分木：ライン実行[%d]行目：x[%d]y[%d] ", i, KData[i].x, KData[i].y);  //確認用テンプレ
    }
  }
  // fprintf(stderr, ":完了[%d]行：x[%d]y[%d]\n",i,KData[i-1].x,KData[i-1].y);
  DO_WHEN_DEBUG(fprintf(stderr, "6\n"));
  cvFree(&pts);
}

/*****************[CV表示：4分木:信頼度2]**************************/
void DRAW_TREE_4c2( IplImage *imgA, CvPoint *pts , int x, int y, int dim) {
  //点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。
  DO_WHEN_DEBUG(fprintf(stderr, "DRAW_TREE_4c2\n"));
  pts[2].x = x;
  pts[2].y = y;
  int radius = 2 * (pts[1].x - pts[0].x);
  if (dim == 0) {
    if (DD == 0) {
      cvCircle(imgA, pts[2], radius, CV_RGB(250, 130, 130), -1, 1, 0);
    }
    else if (DD == 1) {
      cvCircle(imgA, pts[2], radius, CV_RGB(130, 250, 130), -1, 1, 0);
    }
    else if (DD == 2) {
      cvCircle(imgA, pts[2], radius, CV_RGB(130, 130, 250), -1, 1, 0);
    }
    return;
  }

  CvPoint *middle = (CvPoint *) cvAlloc (sizeof (CvPoint) * 10);
  middle->x = pts[0].x + int( ((pts[1].x - pts[0].x) / 2.0) + 0.5 );
  middle->y = pts[0].y + int( ((pts[1].y - pts[0].y) / 2.0) + 0.5 );

  if (middle->x <= x) {
    pts[0].x = middle->x;
  }
  else {
    pts[1].x = middle->x;
  }

  if (middle->y <= y) {
    pts[0].y = middle->y;
  }
  else {
    pts[1].y = middle->y;
  }

  cvFree(&middle);

  DRAW_TREE_4c2( imgA, pts , x, y, dim - 1);
}

/*****************[CV表示：4分木:信頼度1]**************************/
void DRAW_TREE_4c1( IplImage *imgA, CvPoint *pts , int x, int y, int dim) { //点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。
  pts[2].x = x;
  pts[2].y = y;
  int radius = pts[1].x - pts[0].x;
  if (dim == 0) {
    if (DD == 0) {
      cvCircle(imgA, pts[2], radius, CV_RGB(250, 90, 90), -1, 1, 0);
    }
    else if (DD == 1) {
      cvCircle(imgA, pts[2], radius, CV_RGB(90, 250, 90), -1, 1, 0);
    }
    else if (DD == 2) {
      cvCircle(imgA, pts[2], radius, CV_RGB(90, 90, 250), -1, 1, 0);
    }
    return;
  }

  CvPoint *middle = (CvPoint *) cvAlloc (sizeof (CvPoint) * 10);
  middle->x = pts[0].x + int( ((pts[1].x - pts[0].x) / 2.0) + 0.5 );
  middle->y = pts[0].y + int( ((pts[1].y - pts[0].y) / 2.0) + 0.5 );

  if ( middle->x <= x ) {
    pts[0].x = middle->x;
  }
  else {
    pts[1].x = middle->x;
  }

  if (middle->y <= y) {
    pts[0].y = middle->y;
  }
  else {
    pts[1].y = middle->y;
  }

  cvFree(&middle);

  DRAW_TREE_4c1( imgA, pts , x, y, dim - 1);
}


/*****************[CV表示：4分木:ピクセル]**************************/
void DRAW_TREE_4v( IplImage *imgA, CvPoint *pts , int x, int y, int dim) { //点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。
  if (dim == 0) {
    if (DD == 0) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(255, 45, 45), CV_FILLED);
    }
    else if (DD == 1) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(45, 255, 45), CV_FILLED);
    }
    else if (DD == 2) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(45, 45, 255), CV_FILLED);
    }
    return;
  }

  CvPoint *middle = (CvPoint *) cvAlloc (sizeof (CvPoint) * 100);
  middle->x = pts[0].x + int( ((pts[1].x - pts[0].x) / 2.0) + 0.5 );
  middle->y = pts[0].y + int( ((pts[1].y - pts[0].y) / 2.0) + 0.5 );

  if ( middle->x <= x ) {
    pts[0].x = middle->x;
  }
  else {
    pts[1].x = middle->x;
  }

  if (middle->y <= y) {
    pts[0].y = middle->y;
  }
  else {
    pts[1].y = middle->y;
  }

  cvFree(&middle);

  DRAW_TREE_4v( imgA, pts , x, y, dim - 1);
}

/*****************[CV表示：4分木:最深ピクセル]**************************/
void DRAW_TREE_4p( IplImage *imgA, CvPoint *pts , int x, int y, int dim) { //点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。
  if (dim <= 0) {
    if (DD == 0) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(255, 0, 0), CV_FILLED);
    }
    else if (DD == 1) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(0, 255, 0), CV_FILLED);
    }
    else if (DD == 2) {
      cvRectangle(imgA, pts[0], pts[1], CV_RGB(0, 0, 255), CV_FILLED);
    }
    return;
  }

  CvPoint *middle = (CvPoint *) cvAlloc (sizeof (CvPoint) * 100);
  middle->x = pts[0].x + int( ((pts[1].x - pts[0].x) / 2.0) + 0.5 );
  middle->y = pts[0].y + int( ((pts[1].y - pts[0].y) / 2.0) + 0.5 );

  if ( middle->x <= x ) {
    pts[0].x = middle->x;
  }
  else {
    pts[1].x = middle->x;
  }

  if (middle->y <= y) {
    pts[0].y = middle->y;
  }
  else {
    pts[1].y = middle->y;
  }

  cvFree(&middle);

  DRAW_TREE_4p( imgA, pts , x, y, dim - 1);
}

/*****************[CV表示：4分木：ライン]**************************/
void DRAW_TREE_4l( IplImage *imgA, CvPoint *pts , int x, int y, int dim) { //点(線)描画：初期のptsには00と高幅、intｘｙには点座標が入っている。

  if (dim <= 0) {
    return;
  }

  CvPoint *middle = (CvPoint *) cvAlloc (sizeof (CvPoint) * 100);
  middle->x = pts[0].x + int( ((pts[1].x - pts[0].x) / 2.0) + 0.5 );
  middle->y = pts[0].y + int( ((pts[1].y - pts[0].y) / 2.0) + 0.5 );

  if (dim > 0) {
    cvLine( imgA, cvPoint( middle->x, pts[0].y ), cvPoint( middle->x , pts[1].y ),  CV_RGB(0, 0, 0) );
    cvLine( imgA, cvPoint( pts[0].x, middle->y ), cvPoint( pts[1].x, middle->y ),   CV_RGB(0, 0, 0) );
  }
  if ( middle->x <= x ) {
    pts[0].x = middle->x;
  }
  else {
    pts[1].x = middle->x;
  }

  if (middle->y <= y) {
    pts[0].y = middle->y;
  }
  else {
    pts[1].y = middle->y;
  }

  cvFree(&middle);

  DRAW_TREE_4l( imgA, pts , x, y, dim - 1);
}


/****[CV既知点表示：キーボードから直接呼ばれる]****/
void KnownPoint( void *imgA, CvPoint *KData ) { //既知点表示。
  fprintf(stderr, "CVcnt(CV用断面の点群行数)=%d \n", CVcnt); //確認用テンプレ
  IplImage *clone = cvCloneImage( (IplImage *)imgA );

  CvPoint *pts = (CvPoint *) cvAlloc (sizeof (CvPoint) * 200);

  DRAW_TREE_4(clone, KData);//統括にまとめて送って処理してもらう

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  DO_WHEN_DEBUG(fprintf(stderr, "About to show image... "));
  cvShowImage(cv_title, clone);//描画終了後、まとまってから表示。軽量化。
  DO_WHEN_DEBUG(fprintf(stderr, "done.\n"));
  fprintf(stderr, "\nOpenCV再表示完了\n");
  cvReleaseImage(&clone);
  cvFree(&pts);
}

/****[CV 二次元用正規化処理 & 格納]****/
CvPoint *CVset(void) { //data_3dのx/y/zをグローバルなCvPoint型のCVPtsに入れる。
  double x, y;
  if (DD == 0) {
    DP = XX + DX / 2, DM = XX - DX / 2;  //DD=0:X軸断面
  }
  else if (DD == 1) {
    DP = YY + DX / 2, DM = YY - DX / 2;  //DD=1:Y軸断面
  }
  else if (DD == 2) {
    DP = ZZ + DX / 2, DM = ZZ - DX / 2;  //DD=2:Z軸断面
  }
  CvPoint *CVPtsData;
  CVcnt = 0;
  int nuru = 0;
  fprintf(stderr, "\n　確認DD:%d　\n", DD);

  CVPtsData = (CvPoint *) cvAlloc (sizeof (CvPoint) * data_3d.size() / 3); //
  if (CVPtsData == nil) {
    fprintf(stderr, "can't CALLOC CVPtsData.\n");
  }
  for (int i = 0; i < data_3d.size(); i += 3) {
    if (data_3d[i] != nil && data_3d[i + 1] != nil && data_3d[1 + 2] != nil) {  //3次元のXYZがnilでないとき
      if (DM <= data_3d[i + DD] && data_3d[i + DD] < DP) {
        x = ((data_3d[i + DDX] - vox_min[DDX]) / vox_max_divide); //正規格子実験
        y = ((data_3d[i + DDY] - vox_min[DDY]) / vox_max_divide); //正規格子実験
        CVPtsData[CVcnt].x = (int)(x * cv_width);
        CVPtsData[CVcnt].y = cv_height - (int)(y * cv_height); //CVの原点は左上の為。
        if (CVcnt % 100 == 0) {
          fprintf(stderr, "CV%d行GL結果: (%f, %f) \n", CVcnt, data_3d[i + DDX], data_3d[i + DDY] );
        }
        if (CVcnt % 100 == 0) {
          fprintf(stderr, "CV%d行結果: (%d, %d) \n", CVcnt, CVPtsData[CVcnt].x, CVPtsData[CVcnt].y );
        }
        CVcnt++;
      }
    }
    else {
      nuru++;
      if (nuru % 100 == 0) {
        fprintf(stderr, "data_3d=nil %d回目 \r", nuru);
      }
    }
  }
  fprintf(stderr, "\n　結果：%d行読み取り\n", CVcnt);
  CVPts = CVPtsData;
  return(CVPtsData);//cではCVPtsに入る。
  cvFree(&CVPtsData);

}

/*-------------------------------
 [マウス・キーボード・GLUI関係]
 --------------------------------*/

/*****************[マウス(主にCVの画面で使用)]**************************/
void MOUSE(int event, int x, int y, int flags, void *imgA) {
  DO_WHEN_DEBUG(fprintf(stderr, "MOUSE\n"));
  IplImage *clone = cvCloneImage( (IplImage *)imgA );
  CvPoint *pts = (CvPoint *) cvAlloc(sizeof (CvPoint) * 20000); //四分木用のpts領域確保
  static bool MOUSE_FLAG = false;
  int loop_count = 0;
  switch (event) {
    case CV_EVENT_LBUTTONDOWN://左：押しむとフラグが立つ。下記イベントで使用。
      MOUSE_FLAG = true;
      break;
    case CV_EVENT_LBUTTONUP://左：離すとフラグが折れる。
      MOUSE_FLAG = false;
      if (ReviewFlag == 0) {
        DRAW_TREE_4(clone, CVPts);//統括にまとめて送って処理してもらう
        cvShowImage(cv_title, clone);//描画終了後、まとまってから表示。軽量化
        fprintf(stderr, "\nOpenCV再表示完了\n");
      }
      break;
    case CV_EVENT_MBUTTONDOWN://中央：ボクセル分割を再度行う。
      fprintf(stderr, "ボクセル削除開始\n");//確認用
      deleteTVoxcel(root);
      fprintf(stderr, "削除完了\n 再構築開始\n");//確認用
      Make_Voxcel();
      glutPostRedisplay();
      break;
    case CV_EVENT_RBUTTONDOWN://右:一点ずつ配置、もしくは削除
      if (DeleteFlag == 0) { //DeleteFlagが０の時、追加
        point_plus((int)x, (int)y, imgA);
        if (ReviewFlag == 0) {
          DRAW_TREE_4(clone, CVPts);//統括にまとめて送って処理してもらう
          cvShowImage(cv_title, clone);//描画終了後、まとまってから表示。軽量化
          fprintf(stderr, "\nOpenCV再表示完了\n");
        }
      }
      else if (DeleteFlag == 1) {   //DeleteFlag==1の時、削除
        double DGLp, Dxyz;//Deleate_GL_Point(各軸の断面中央値),Deleate_XYZ(各軸の断面幅)
        DGLp = 0.0;
        Dxyz = 0.0;
        if (DD == 0) {
          DGLp = XX;
          Dxyz = DX; //DD=0:X軸断面
        }
        else if (DD == 1) {
          DGLp = YY;
          Dxyz = DX; //DD=1:Y軸断面
        }
        else if (DD == 2) {
          DGLp = ZZ;
          Dxyz = DX; //DD=2:Z軸断面
        }
        double DqtreeX = cv_width / (pow(2, (double)vox_value)); //CVの削除する幅:DimensionQuadtreeX
        double DqtreeY = cv_height / (pow(2, (double)vox_value)); //CVの削除する幅:DimensionQuadtreeY.
        double DqtX_count, DqtY_count; //CVの削除するピクセルの下限:QuadtreeX,Y
        for (DqtX_count = 0; DqtX_count < x; DqtX_count += DqtreeX)
        {}
        DqtX_count -= DqtreeX; //削除するX軸左辺を特定
        for (DqtY_count = 0; DqtY_count < y; DqtY_count += DqtreeY)
        {}
        DqtY_count -= DqtreeY; //削除するY軸上辺を特定
        //CVの点・窓を削除、再描画
        for (int i = 0; i < CVcnt; i++) {
          if (DqtX_count <= CVPts[i].x && CVPts[i].x < DqtX_count + DqtreeX) {
            if (DqtY_count <= CVPts[i].y && CVPts[i].y < DqtY_count + DqtreeY) {
              loop_count++;
              CVPts[i].x = nil;
              CVPts[i].y = nil;
              if (i % 10 == 0) {
                fprintf(stderr, "CVPts[%d]:[%d,%d]ライン削除完了\n", loop_count, x, y); //確認用
              }
            }
          }
        }
        /*for(int i=0; i<CVcnt; i++){//CV画像書きなおし
         pts[0] = cvPoint(0,0);
         pts[1] = cvPoint(clone->width,clone->height);
         //if(CVPts[i].x != nil)DRAW_TREE_4( clone, pts, CVPts[i].x, CVPts[i].y, dimension);
         //if(CVPts[i].x != nil)DRAW_TREE_4( clone, pts, CVPts[i].x, CVPts[i].y, vox_value);
         }*/
        if (ReviewFlag == 0) {
          DRAW_TREE_4(clone, CVPts);//統括にまとめて送って処理してもらう
          cvShowImage(cv_title, clone);
          fprintf(stderr, "\nOpenCV再表示完了\n");
        }

        //GLの点を削除
        double nearX = vox_max_divide * (((double)x) / cv_width) + vox_min[DDX]; //指定点のGL側:CVのX座標相当
        double nearY = vox_max_divide * (((double)(cv_height - y)) / cv_height) + vox_min[DDY]; //指定点のGL側：CVのY座標
        double Dotree = vox_max_divide / (pow(2, (double)vox_value)); //GLの削除するボクセルの幅
        double DotX_count, DotY_count; //GLの削除するボクセルの下限（CVのXとYに相当）
        for (DotX_count = vox_min[DDX]; DotX_count < nearX; DotX_count += Dotree)
        {}
        DotX_count -= Dotree; //削除するX座標左辺を特定
        for (DotY_count = vox_min[DDY]; DotY_count < nearY; DotY_count += Dotree)
        {}
        DotY_count -= Dotree; //削除するY座標上辺を特定

        for (int i = 0; i < data_3d.size(); i += 3) {
          if ((DGLp - Dxyz / 2) <= data_3d[i + DD] && data_3d[i + DD] < (DGLp + Dxyz / 2)) { //採取断面の幅内において
            if (DotX_count <= data_3d[i + DDX] && data_3d[i + DDX] < DotX_count + Dotree) {
              if (DotY_count <= data_3d[i + DDY] && data_3d[i + DDY] < DotY_count + Dotree) {
                fprintf(stderr, "第%d行目削除,[%f,%f,%f]\n", i / 3, data_3d[i + DD], data_3d[i + DDX], data_3d[i + DDY]); //確認用
                data_3d[i] = data_3d[i + 1] = data_3d[i + 2] = nil;   //xyz座標全てにnilを代入
              }
            }
          }
        }
        glutPostRedisplay();
        fprintf(stderr, "Rクリック：処理完了 \n");//最終表示
      }
      else {
        fprintf(stderr, "異常事態です。DeleteFlagを確認して下さい \n");//異常
      }
      break;
    default:
      break;
  }

  if (event == CV_EVENT_MOUSEMOVE && MOUSE_FLAG == true) {
    if (DeleteFlag == 0) { //DeleteFlagが０の時、追加
      point_plus((int)x, (int)y, imgA);
    }
  }
  cvReleaseImage(&clone);
  cvFree(&pts);
}

/*****************[キーボード(主にGLの窓で使用)]**************************/

void keyboard(unsigned char key, int x, int y) {
  switch (key) { //[a, c, d, e, f, g, h, i, j, k, l, m, n, o, p, r, w, x, z,Esc]が使用中
    case ';':
      fprintf(stderr, "key ; was hit.");
      break;
      //    [GL用]
    case 'X'://X軸回転
      rotX -= 2;
      glutPostRedisplay();
      break;
    case 'x'://X軸回転
      rotX += 2;
      glutPostRedisplay();
      break;
    case 'Y'://Y軸回転
      rotY -= 2;
      glutPostRedisplay();
      break;
    case 'y'://Y軸回転
      rotY += 2;
      glutPostRedisplay();
      break;
    case 'Z'://Z軸回転
      rotZ -= 2;
      glutPostRedisplay();
      break;
    case 'z'://Z軸回転
      rotZ += 2;
      glutPostRedisplay();
      break;
    case 'f'://カメラ後退
      dist += 1;
      glutPostRedisplay();
      break;
    case 'F'://カメラ前進
      dist -= 1;
      glutPostRedisplay();
      break;
    case 'd'://GLデフォルト表示
      dist = 100;
      rotX = rotY = rotZ = 0;
      fprintf(stderr, "XYZデフォルト位置\n"); //確認用
      glutPostRedisplay();
      break;
    case 'D'://GL表示
      glutPostRedisplay();
      fprintf(stderr, "画像再表示\n"); //確認用
      break;
    case 'e'://採取断面軸変更
      DD++;
      if (DD > 2) {
        DD = 0;
      }
      if (DD == 0) {
        DDX = 1, DDY = 2;  //DD=0:X軸断面
      }
      else if (DD == 1) {
        DDX = 2, DDY = 0;  //DD=1:Y軸断面
      }
      else if (DD == 2) {
        DDX = 0, DDY = 1;  //DD=2:Z軸断面
      }
      fprintf(stderr, "DD:%d (X:0, Y:1, Z:2)\n", DD); //確認用
      break;
    case 'E'://採取断面軸変更
      DD--;
      if (DD < 0) {
        DD = 2;
      }
      if (DD == 0) {
        DDX = 1, DDY = 2;  //DD=0:X軸断面
      }
      else if (DD == 1) {
        DDX = 2, DDY = 0;  //DD=1:Y軸断面
      }
      else if (DD == 2) {
        DDX = 0, DDY = 1;  //DD=2:Z軸断面
      }
      fprintf(stderr, "DD:%d (X:0, Y:1, Z:2)\n", DD); //確認用
      break;
    case 'a'://表示ボクセル階層数減少
      if (vox_value > 0) {
        vox_value -= 1;
        DAxis(0);//断面幅変更
        if (dVFlag == 1) {
          value_counter_con();  //点群信頼度のMAX変算出
        }
        fprintf(stderr, "表示階層/最深階層[%d/%d]:断面幅[%f]\n", vox_value, MAX_DEVIDE_FREQUENCY, DX);
        fprintf(stderr, "採取軸確認、X：%d、Y:%d、Z:%d：\n", countDX, countDY, countDZ);
        glutPostRedisplay();
      }
      else {
        fprintf(stderr, "ボクセル階層変更不可：深⇒浅：最浅階層表示中\n");
      }
      break;
    case 'A'://表示ボクセル階層数増加
      if (vox_value < MAX_DEVIDE_FREQUENCY) {
        vox_value += 1; //ボクセル階層変更
        DAxis(1);//断面幅変更
        if (dVFlag == 1) {
          value_counter_con();  //点群信頼度のMAX変算出
        }
        fprintf(stderr, "表示階層/最深階層[%d/%d]:断面幅[%f]\n", vox_value, MAX_DEVIDE_FREQUENCY, DX);
        fprintf(stderr, "採取軸確認、X：%d、Y:%d、Z:%d：\n", countDX, countDY, countDZ);
        glutPostRedisplay();
      }
      else {
        fprintf(stderr, "ボクセル階層変更不可：浅⇒深：最深階層表示中\n");
      }
      break;
    case 'g'://ボクセル描画方法変更：ワイヤー&採集軸
      dVFlag = 0;
      fprintf(stderr, "ボクセル表示方法変更：基本 \n");//確認
      break;
    case 'G'://ボクセル描画方法変更：信頼度
      dVFlag = 1;
      value_counter_con();//指定階層のボクセルの内、内包点群最大の個数計算
      fprintf(stderr, "ボクセル表示方法変更：信頼度 \n"); //確認
      break;
    case 'l'://
      acFlag = 0;
      fprintf(stderr, "GL：acFlagOFF：全ボクセルを表示：\n");
      break;
    case 'L'://
      acFlag = 1;
      fprintf(stderr, "GL：acFlagON：任意断面のみの表示:：\n");
      break;
      //    [CV用：主に断面関係]
    case 'i'://採取X軸移動
      XX -= DX;
      countDX -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取X軸移動、断面幅の小側現在値：%d：\n", countDX);
      break;
    case 'I'://採取X軸移動
      XX += DX;
      countDX += pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取X軸移動、断面幅の小側現在値：%d：\n", countDX);
      break;
    case 'j'://採取y軸移動
      YY -= DX;
      countDY -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取Y軸移動、断面幅の小側現在値：%d：\n", countDY);
      break;
    case 'J'://採取Y軸移動
      YY += DX;
      countDY += pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取Y軸移動、断面幅の小側現在値：%d：\n", countDY);
      break;
    case 'k'://採取z軸移動
      ZZ -= DX;
      countDZ -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取Z軸移動、断面幅の小側現在値：%d：\n", countDZ);
      break;
    case 'K'://採取Z軸移動
      ZZ += DX;
      countDZ += pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value));
      glutPostRedisplay();
      fprintf(stderr, "採取Z軸移動、断面幅の小側現在値：%d：\n", countDZ);
      break;
    case 'c'://断面のx・yに当たる値を取得
      DP = XX + DX / 2;
      DM = XX - DX / 2;
      CVPts = CVset();
      fprintf(stderr, "cv用データ入力完了 \n");
      break;
    case 'C'://OpenCVの窓に断面の点群を生成
      fprintf(stderr, "cv描画開始 \n");
      KnownPoint( imgA, CVPts );//既知点の四分木&表示。cで作られるファイルに対応＝データ型を使用
      fprintf(stderr, "cv描画完了 \n");
      break;
    case 'm'://マウス操作時、点を追加
      DeleteFlag = 0;
      fprintf(stderr, "DeleteFlag変更：点追加待機\n");
      break;
    case 'M'://マウス操作時、点を削除
      DeleteFlag = 1;
      fprintf(stderr, "DeleteFlag変更：点削除待機\n");
      break;
    case 'n'://点を追加時、断面中央
      PlusFlag = 0;
      fprintf(stderr, "PlusFlag変更：点追加時⇒1点ずつ(0:0/1)\n");
      break;
    case 'N'://点を追加時、全最深ボクセル
      PlusFlag = 1;
      fprintf(stderr, "PlusFlag変更：点追加時⇒対応全最深ボクセル(1:0/1)\n");
      break;
    case 'o'://点群表示する際、ピクセル表示するか⇒する
      PixelFlag = 0;
      fprintf(stderr, "PixelFlag変更：CV表示時、ボクセル着色『有』\n");
      break;
    case 'O'://点群表示する際、ピクセル表示するか⇒しない
      PixelFlag = 1;
      fprintf(stderr, "PixelFlag変更：CV表示時、ボクセル着色『無』\n");
      break;
    case 'p'://点群表示する際、信頼度を表示するか⇒する
      ConFlag = 0;
      fprintf(stderr, "PixelFlag変更：CV表示時、信頼度段階表示『有』\n");
      break;
    case 'P'://点群表示する際、信頼度を表示するか⇒しない
      ConFlag = 1;
      fprintf(stderr, "PixelFlag変更：CV表示時、信頼度段階表示『無』\n");
      break;
    case 'q'://マウス操作を行った際、再表示を行うか⇒行う
      ReviewFlag = 0;
      fprintf(stderr, "PixelFlag変更：CVマウス時、再表示『有』\n");
      break;
    case 'Q'://マウス操作を行った際、再表示を行うか⇒行わない
      ReviewFlag = 1;
      fprintf(stderr, "PixelFlag変更：CVマウス時、再表示『無』\n");
      break;
      //    [ファイル入出力用]
    case 'w'://ファイル書き込み
      FileOutput();
      break;
    case 'W'://ファイル書き込み予定のデータ(data_3d)の現在値閲覧
      DataRead();
      break;
    case 'r'://ファイル読み取り:
      ReFileInput(OUTPUTDATA);
      break;
    case 'R'://ファイル読み取り：元データ
      ReFileInput(SCAN3D_POINT_DATA);
      break;
    case 'h'://ヘルプ
      fprintf(stderr, "***キーボードヘルプ***\n");
      fprintf(stderr, " a/A:GL:階層浅/深,c:断面取得,C:断面表示\n d:GL再表示,e/E:採取断面変更+/-,f/F:ｶﾒﾗ後/前\n");
      fprintf(stderr, " g/G:GLのボクセルカラー:基本/信頼度\n h:ｷｰﾍﾙﾌﾟ,H:その他確認,I/J/K:断面X/Y/Z軸移動\n");
      fprintf(stderr, " l/L:GLのボクセル全体or任意断面\n m/M:CV:追加/削除,n/N:CV:点追加時,断面中央/全最深ボクセルの同X/Y\n");
      fprintf(stderr, " 0/O:CV表示:ピクセルON/OFF,p/P:CV表示：信頼度ON/OFF,q/Q:CVマウス時再表示：ON/OF\n");
      fprintf(stderr, " r:出力ﾌｧｲﾙ読取,R:元ﾌｧｲﾙ読取,\n w:ﾌｧｲﾙ書込,W:ﾌｧｲﾙ書込内容閲覧,x/y/z:カメラ移動\n");
      fprintf(stderr, " Esc:デバッグ終了\n");
      break;
    case 'H'://その他確認
      fprintf(stderr, "***その他のヘルプ***\n");
      fprintf(stderr, " 色相信頼度：RED<YELLOW<GREEN<BULE\n");
      fprintf(stderr, " 表示階層/最深階層[%d/%d]:断面幅[%f]\n", vox_value, MAX_DEVIDE_FREQUENCY, DX);
      fprintf(stderr, " GL:採取軸確認:%d (X:0, Y:1, Z:2)\n", DD);
      fprintf(stderr, " GL:採取軸位置確認、X：%d、Y:%d、Z:%d：\n", countDX, countDY, countDZ);
      fprintf(stderr, " CV:点追加/削除:%d,点追加1点/全最深ボクセル:%d,\n", DeleteFlag, PlusFlag);
      break;
    case 27:
      exit(0);
      break;
    default:
      break;

  }
}

/*****************[3次元データ現在値]**************************/
void DataRead(void) {
  int loopcnt = 0;
  int nullcnt = 0;
  for (int i = 0; i < data_3d.size(); i += 3) {
    if (data_3d[i] == nil) {
      if (nullcnt % 100 == 0) {
        fprintf(stderr, "data_3d[%d]行目データ=nil：%f, %f,%f \n", loopcnt, data_3d[i], data_3d[i + 1], data_3d[i + 2]);
      }
      nullcnt++;
    }
    if (loopcnt % 10000 == 0) {
      fprintf(stderr, "data_3d[%d]行目データ：%f, %f,%f \n", loopcnt, data_3d[i], data_3d[i + 1], data_3d[i + 2]);
    }
    loopcnt++;
  }
  fprintf(stderr, "data_3dに存在する点群の行数：[%d]\n", loopcnt - nullcnt);
}


/*************************
 [点群信頼度、CONから追加]
 ***************************/

/*****************[ボクセル内点群カウント指示]********************/
void value_counter_con(void) {
  confidence_max = 0;
  counterCpoint(root, vox_value);
  fprintf(stderr, "confidence_max:%d\n", confidence_max); //前のものに上書きするためには\r
}


/*****************[ボクセル内点群カウント]********************/
void counterCpoint(TVoxcel *voxcel, int value) {
  if (voxcel->value == value) {
    confidence_p = 0;
    if (confidence_max < voxcel->point_number) {
      confidence_max = voxcel->point_number;  //より多くの内包点になれば、そちらに変更
    }
  }
  else {   //再帰：ボクセルのvalueがvalue（指定した階層）ではなく、さらに〜分木の構成の場合、
    if (voxcel->child != nil) { //N分木ならN回、自分の子を行う。例えば5階層のボクセル一覧を表示するなら、5階層までもぐって実行
      for (int i = 0; i < voxcel->child_size; i++) {
        counterCpoint(voxcel->child[i], value);
      }
    }
  }
}

/*************************
 [ファイル入出力]
 ***************************/


/*****************[ファイル出力]**************************/
void FileOutput(void) {
  int FOp_count = 0;
  fprintf(stderr, "書き込み処理->\n");
  FILE *fp = fopen(OUTPUTDATA, "w");
  fprintf(fp, "#TVOXCEL:level%d\n", MAX_DEVIDE_FREQUENCY);  //stderrではなく、ファイル一行目に書き込み
  if (fp == nil) {                     // オープンに失敗した場合
    fprintf(stderr, "cannot open:OutputFILE\n");         // エラーメッセージを出して
    exit(1);                       // 異常終了
  }
  for (int i = 0; i < data_3d.size(); i += 3) {           //inputファイルの行数回数繰り返す
    if (data_3d[i] != nil && data_3d[i + 1] != nil && data_3d[i + 2] != nil) { //入っているデータが全てnilではない＝値が入っている時
      if (FOp_count % 10000 == 0) {
        fprintf(stderr, "output[%d]行目データ：%f, %f,%f \r", FOp_count, data_3d[i], data_3d[i + 1], data_3d[i + 2]);  //確認用
      }
      fprintf(fp, "%f %f %f \n", data_3d[i], data_3d[i + 1], data_3d[i + 2]); //
      FOp_count++;
    }
  }
  fprintf(stderr, "\n [%d]行のデータ書き込み", FOp_count);
  fprintf(stderr, "\n >書き込み完了\n");

  fclose(fp);
}

/*****************[ファイル再入力]**************************/

void ReFileInput(const char *Input_data) { //その場でアウトプットしたファイルを読み込む為の物

  for (int i = 0; i < data_3d.size(); i += 3) { //inputファイルの行数回数繰り返す
    data_3d[i] = data_3d[i + 1] = data_3d[i + 2] = nil; //全点にnil
  }
  int RFI_count = 0;
  int size_count = 0;
  fprintf(stderr, "点データの読み込み開始\n");
  FILE *fp = fopen(Input_data, "r");
  double tx, ty, tz;
  char tmp_c[255];
  while (fgets(tmp_c, 100, fp) && !feof(fp)) {
    if (tmp_c[0] != 0x23) {   //先頭行が#でなければ読み込む
      sscanf(tmp_c, "%lf %lf %lf\n", &tx, &ty, &tz);
      if (data_3d.size() == 0 || vox_min[0] > tx) {
        vox_min[0] = tx;
      }
      if (data_3d.size() == 0 || vox_max[0] < tx) {
        vox_max[0] = tx;
      }
      if (data_3d.size() == 0 || vox_min[1] > ty) {
        vox_min[1] = ty;
      }
      if (data_3d.size() == 0 || vox_max[1] < ty) {
        vox_max[1] = ty;
      }
      if (data_3d.size() == 0 || vox_min[2] > tz) {
        vox_min[2] = tz;
      }
      if (data_3d.size() == 0 || vox_max[2] < tz) {
        vox_max[2] = tz;
      }
      if (data_3d[size_count] != nil || data_3d[size_count + 1] != nil || data_3d[size_count + 2] != nil) {
        data_3d.push_back(tx), data_3d.push_back(ty), data_3d.push_back(tz);//いずれかに値が入っていれば新たに『追加』する
      }
      else {   //いずれもnilの場合（前回のdata_3dが消されている場合、上書きする
        data_3d[size_count] = tx;
        data_3d[size_count + 1] = ty;
        data_3d[size_count + 2] = tz;
      }

      if (RFI_count % 100000 == 0)fprintf(stderr, "%d行目まで新規読み込み(%f, %f, %f)\r",
                                            RFI_count, data_3d[size_count], data_3d[size_count + 1], data_3d[size_count + 2]);
      RFI_count++;
      size_count += 3;
    }
  }
  //対象オブジェクトの最大幅を決定
  if (vox_max[0] - vox_min[0] >= vox_max[1] - vox_min[1] && vox_max[0] - vox_min[0] >= vox_max[0] - vox_min[0]) {
    vox_max_divide = vox_max[0] - vox_min[0];
  }
  else if (vox_max[1] - vox_min[1] >= vox_max[2] - vox_min[2]) {
    vox_max_divide = vox_max[1] - vox_min[1];
  }
  else {
    vox_max_divide = vox_max[2] - vox_min[2];
  }
  fclose(fp);
  fprintf(stderr, "\n [%d]点のデータを読み込み完了\n", RFI_count);
  fprintf(stderr, "ボクセル削除開始\n");//確認用
  deleteTVoxcel(root);
  fprintf(stderr, "削除完了\n 再構築開始\n");//確認用
  Make_Voxcel();
  glutPostRedisplay();
}


void read_Point(const char *Input_data) {
  int rP_count = 0;
  std::cerr << "Start reading pointcloud." << std::endl;
  FILE *fp = fopen(Input_data, "r");
  double tx, ty, tz;
  char tmp_c[255];
  while (fgets(tmp_c, 100, fp) && !feof(fp)) { // ERROR: DO NOT USE FEOF.
    if (tmp_c[0] != '#') {
      sscanf(tmp_c, "%lf %lf %lf\n", &tx, &ty, &tz);  // WARNING: SHOULD USE std::cin.
      // WARNING: THE FOLLOWING CODE SHOULD BE REPLACED WITH MAX.
      if (data_3d.size() == 0 || vox_min[0] > tx) {
        vox_min[0] = tx;
      }
      if (data_3d.size() == 0 || vox_max[0] < tx) {
        vox_max[0] = tx;
      }
      if (data_3d.size() == 0 || vox_min[1] > ty) {
        vox_min[1] = ty;
      }
      if (data_3d.size() == 0 || vox_max[1] < ty) {
        vox_max[1] = ty;
      }
      if (data_3d.size() == 0 || vox_min[2] > tz) {
        vox_min[2] = tz;
      }
      if (data_3d.size() == 0 || vox_max[2] < tz) {
        vox_max[2] = tz;
      }
      data_3d.push_back(tx);
      data_3d.push_back(ty);
      data_3d.push_back(tz);
      if (rP_count % 100000 == 0) {
        std::cerr << "Read " << rP_count << std::endl;
      }
      rP_count++;
    }
  }
  if (vox_max[0] - vox_min[0] >= vox_max[1] - vox_min[1]
      && vox_max[0] - vox_min[0] >= vox_max[0] - vox_min[0]) {
    vox_max_divide = vox_max[0] - vox_min[0];
  }
  else if (vox_max[1] - vox_min[1] >= vox_max[2] - vox_min[2]) {
    vox_max_divide = vox_max[1] - vox_min[1];
  }
  else {
    vox_max_divide = vox_max[2] - vox_min[2];
  }
  fclose(fp);
  std::cerr << "Finished reading pointcloud." << std::endl;
}

/*****************[ボクセル作成]**************************/

void Make_Voxcel(void) {
  //ルートボクセル初期化
  int MV_count = 0;
  fprintf(stderr, "\nボクセルデータ作成開始\n");
  root = new TVoxcel();//root＝根＝xyz全ての最大最小が入った大本のボックスセル
  initTVoxcel(root);
  root->min_x = vox_min[0];
  root->min_y = vox_min[1];
  root->min_z = vox_min[2];
  root->max_x = vox_min[0] + vox_max_divide;
  root->max_y = vox_min[1] + vox_max_divide;
  root->max_z = vox_min[2] + vox_max_divide;
  for (int i = 0; i < data_3d.size(); i += 3) {
    if (data_3d[i] != nil && data_3d[i + 1] != nil && data_3d[1 + 2] != nil) {
      devideTVoxcelByPoint(root, data_3d[i], data_3d[i + 1], data_3d[i + 2], MAX_DEVIDE_FREQUENCY, 2);
      if (MV_count % 100000 == 0) {
        fprintf(stderr, "%d行目までのボクセル構築完了\r", MV_count);
      }
      MV_count++;
    }
  }
  fprintf(stderr, "\nボクセルデータ作成完了(%d)\n", MV_count);
}


/**-----------------------------
 [その他]
 -------------------------------**/

/*****************[配色する断面幅変更]**************************/
void DAxis(int Flag) {
  if (Flag == 0) { //小ボクセル→大ボクセル
    if (countDX % (int)pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value))  == 0) {
      XX += DX / 2;
    }
    else {
      XX -= DX / 2;
      countDX -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value - 1));
    }
    if (countDY % (int)pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value)) == 0) {
      YY += DX / 2;
    }
    else {
      YY -= DX / 2;
      countDY -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value - 1));
    }
    if (countDZ % (int)pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value)) == 0) {
      ZZ += DX / 2;
    }
    else {
      ZZ -= DX / 2;
      countDZ -= pow(2, double(MAX_DEVIDE_FREQUENCY - vox_value - 1));
    }
  }
  else if (Flag == 1) { //大ボクセル→小ボクセル
    XX -= DX / 4;
    YY -= DX / 4;
    ZZ -= DX / 4;
  }
  DX = vox_max_divide / (pow(2, (double)vox_value));
}

/*****************[CVによる点追加時]**************************/
void point_plus(int plus_x, int plus_y, void *imgA) {
  int CVX[2] = {0, 0}, CVY[2] = {0, 0}; //指定した際のピクセルXY軸の上底下底
  CVPts[CVcnt].x = (int)plus_x; //最新のCVの窓の値が格納されている。
  CVPts[CVcnt].y = (int)plus_y; //最新のCVの窓の値が格納されている。
  //ピクセルの上限下限を特定
  for (int i = 0; i < plus_x; i += pow(2.0, (double)vox_value)) {
    CVX[0] = i;
    CVX[1] = i + pow(2, (double)(MAX_DEVIDE_FREQUENCY - vox_value));
  }
  for (int i = cv_height; i > plus_y; i -= pow(2.0, (double)vox_value)) {
    CVY[0] = cv_height - i;
    CVY[1] = cv_height - (i - pow(2, (double)(MAX_DEVIDE_FREQUENCY - vox_value)));
  }
  fprintf(stderr, "ピクセル各上限下限X{%d,%d},{%d,%d}\n", CVX[0], CVX[1], CVY[0], CVY[1]); //確認用
  //断面中央に一点追加
  if (PlusFlag == 0) {
    //cvtoglpointplus((double)plus_x, (double)plus_y,0);
    if (DD == 0) { //X断面時の格納
      glx = XX; //X断面の中央を格納
      gly = ((double)plus_x / (double)cv_width) * vox_max_divide + vox_min[1]; //CVのｘ軸からGLのY軸に変換
      glz = (double)(cv_height - plus_y) / (double)cv_height * vox_max_divide + vox_min[2]; //CVのcv_height-ｙ軸からGLのZ軸に変換
    }
    else if (DD == 1) { //Y断面時の格納
      glx = ((double)(cv_height - plus_y) / (double)cv_height) * vox_max_divide + vox_min[0]; //CVのcv_height-ｙ軸からGLのZ軸に変換
      gly = YY;
      glz = ((double)plus_x / (double)cv_width) * vox_max_divide + vox_min[2]; //CVのｘ軸からGLのY軸に変換
    }
    else if (DD == 2) { //Z断面時の格納
      glx = ((double)(plus_x) / (double)cv_width) * vox_max_divide + vox_min[0]; //CVのｘ軸からGLのY軸に変換
      gly = ((double)(cv_height - plus_y) / (double)cv_height) * vox_max_divide + vox_min[1]; //CVのcv_height-ｙ軸からGLのZ軸に変換
      glz = ZZ;
    }
    fprintf(stderr, "kakunin： ( %d, %d) \n",  plus_x, plus_y);//確認用
    data_3d.push_back(glx), data_3d.push_back(gly), data_3d.push_back(glz);//点群のファイル型に対し、点の追加
    fprintf(stderr, "CV追加点：%d: (%d, %d) \n", CVcnt, CVPts[CVcnt].x, CVPts[CVcnt].y );//確認用
    fprintf(stderr, "GL追加点： (%f, %f, %f) \n", data_3d[data_3d.size() - 3], data_3d[data_3d.size() - 2], data_3d[data_3d.size() - 1]); //確認用
    CVcnt++;
  }
  //指定ピクセルに相当するボクセル内の全最深ボクセル中央に点を追加
  else if (PlusFlag == 1) {
    if (DD == 0) { //X断面時の格納
      gly = ((double)plus_x / (double)cv_width) * vox_max_divide + vox_min[1]; //CVのｘ軸からGLのY軸に変換
      glz = (double)(cv_height - plus_y) / (double)cv_height * vox_max_divide + vox_min[2]; //CVのcv_height-ｙ軸からGLのZ軸に変換
    }
    else if (DD == 1) { //Y断面時の格納
      glx = ((double)(cv_height - plus_y) / (double)cv_height) * vox_max_divide + vox_min[0]; //CVのcv_height-ｙ軸からGLのZ軸に変換
      glz = ((double)plus_x / (double)cv_width) * vox_max_divide + vox_min[2]; //CVのｘ軸からGLのY軸に変換
    }
    else if (DD == 2) { //Z断面時の格納
      glx = ((double)(plus_x) / (double)cv_width) * vox_max_divide + vox_min[0]; //CVのｘ軸からGLのY軸に変換
      gly = ((double)(cv_height - plus_y) / (double)cv_height) * vox_max_divide + vox_min[1]; //CVのcv_height-ｙ軸からGLのZ軸に変換
    }
    for (double i = DM; i < DP; i += DX / pow(2, (double)MAX_DEVIDE_FREQUENCY)) {
      if (DD == 0) {
        glx = i;
      }
      else if (DD == 1) {
        gly = i;
      }
      else if (DD == 2) {
        glz = i;
      }
      data_3d.push_back(glx), data_3d.push_back(gly), data_3d.push_back(glz);
      fprintf(stderr, "CV追加点：%d: (%d, %d) \n", CVcnt, CVPts[CVcnt].x, CVPts[CVcnt].y );//確認用
      fprintf(stderr, "GL追加点： (%f, %f, %f) \n", data_3d[data_3d.size() - 3], data_3d[data_3d.size() - 2], data_3d[data_3d.size() - 1]); //確認用
      CVcnt++;
    }
  }
  fprintf(stderr, "点群追加処理完了 \n");//最終表示
}




/*---------------
 GLUI実験中
 ---------------*/
/*****************[終了ボタン]**************************/
void gluiCallbackExit(int num) {
  exit(0);
}

/*****************[デフォルトボタンver.01]**************************/
void gluiCallbackDef(int num) {
  gluiXYZ[0] = 0.0;
  gluiXYZ[1] = 0.0;
  gluiXYZ[2] = 0.0;
  for (int i = 0; i < 16; i++) {
    if (i == 0 || i == 5 || i == 10 || i == 15) {
      rotate[i] = 1;
    }
    else {
      rotate[i] = 0;
    }
  }
  gluiXYZ[0] = 0;
  gluiXYZ[1] = 0;
  gluiXYZ[2] = 0;
}
