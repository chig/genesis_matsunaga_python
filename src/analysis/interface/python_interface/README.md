% GENESIS 解析ツール群における Python インターフェース機能追加 関連ドキュメント
% 2025.6

# 概要

* ダウンロード・コンパイル方法
* 回帰テスト実行方法
* ディレクトリ・ファイル構造

# 必要環境

```
gfortran 10（GCC 10） 以降が対応
```

# ダウンロード・コンパイル方法

* github から GENESIS をダウンロードしてから、作業用ブランチ develop に変更

```
$ git clone git@github.com:matsunagalab/genesis.git
$ cd genesis/
$ git checkout develop
```

* Python環境を構築

```
$ brew install uv
$ uv venv --python=python3.11
$ source .venv/bin/activate
$ uv pip install torch torchvision torchaudio nglview numpy mdtraj MDAnalysis plotly jupyterlab py3Dmol scikit-learn
```

* mbar_analysis と msd_analysis の module 名の名前競合を回避するために、src/analysis/interface/mbar_analysis に mbar_analysis 関連の module 名を置換したファイル生成 (automake などの前に必要)

```
$ cd src/analysis/interface/python_interface
$ python mbar_rename.py
```

* GENESIS をインストールすると同時に、Python interface もコンパイルされる (LAPACK は必要)。

```
$ cd genesis
$ autoscan
$ autoheader
$ mkdir m4
$ aclocal
$ autoconf
$ libtoolize
$ automake -a
$ ./configure LAPACK_LIBS="-L/usr/local/lib -llapack -lblas"
$ make
$ make install
# 環境変数
$ cd lib/
$ export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
$ cd ../src/analysis/interface/python_interface/
$ export PYTHONPATH=$(pwd):$PYTHONPATH
$ cd ../../..
```

# Jupyter Notebook 実行

```
$ .venv/bin/jupyterlab
# demo.ipynb を開く
```

# 回帰テスト実行方法

```
$ cd genesis/src/analysis/interface/python_interface
$ ./all_run.sh
```

all_run.sh の内容は以下の通り。全回帰テストを実行する。

```
#!/bin/bash

python crd_convert.py
python trj_analysis.py
python wham_analysis.py
python mbar_analysis_umbrella_1d.py
python mbar_analysis_umbrella_block.py
python avecrd_analysis.py
python kmeans_clustering.py
python hb_analysis_count_atom.py
python hb_analysis_count_snap.py
python rmsd_analysis.py
python drms_analysis.py
python rg_analysis.py
python msd_analysis.py
python diffusion_analysis.py
python test_mdanalysis.py
python test_mdtraj.py
```

または、個別に解析ツールを実行したいときは、例えば

```
$ python rmsd_analysis.py
```

など。

## Python スクリプトの編集

必要に応じて編集する。例えば rmsd_analysis.py の場合
* s_molecule を生成するための PDB/PSF ファイルのパスを書く : pdb_path / psf_path
* s_trajectory を生成するための crd_convert 実行の(inp ファイルの)キーワードを crd_convert の引数に書く
* 解析のための(inp ファイルの)キーワードを trj_analysis の引数に書く

```
import os
import pathlib
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


def test_rmsd_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    with genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = "BPTI_run.dcd",
                    md_step = 10,
                    mdout_period = 1,
                    ana_period = 1,
                    repeat = 1,
                    ),
                ],
            trj_format = "DCD",
            trj_type = "COOR+BOX",
            trj_natom = 0,
            selection_group = ["all", ],
            fitting_method = "NO",
            fitting_atom = 1,
            check_only = False,
            pbc_correct = "NO",
            ) as trajs:
        for t in trajs:
            d = genesis_exe.rmsd_analysis(
                    mol, t,
                    selection_group = ["sid:BPTI and an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    analysis_atom  = 1,
                    )
            print(d.rmsd, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_rmsd_analysis()


if __name__ == "__main__":
    main()
```

* trj_analysis.py スクリプトは、CustomTestCase クラスを使って回帰テストを実装したため crd_convert のキーワードを custom_test_case.py で書いている。

# ディレクトリ・ファイル構造

## ディレクトリ

* genesis/src/analysis/interface/python_interface に全ファイルが存在 (以下で説明)
* genesis/src/analysis/interface/mbar_analysis に mbar_rename.py によって作成される mbar_analysis 関連の module 名変更ファイルが存在 (module名変更しただけなので内容説明省略)

## ファイル

### s_molecule 構造体関連

|ファイル名           |内容                                         |
|:--------------------|:--------------------------------------------|
|s_molecule.py        |Python s_molecule クラス定義、Python s_molecule/s_molecule_c 相互変換、入力ファイルから Python s_molecule オブジェクト生成/解放|
|s_molecule_c.py      |s_molecule_c クラス定義                      |
|s_molecule_c_mod.fpp |s_molecule_c/s_molecule 相互変換、s_molecule_c メモリ確保/解放|
|define_molecule.fpp  |入力ファイルから s_molecule_c 作成           |

### s_trajectories 構造体関連

|ファイル名               |内容                                           |
|:------------------------|:----------------------------------------------|
|s_trajectories.py        |Python s_trajectory クラス定義                 |
|s_trajectories_c.py      |s_trajectory_c クラス定義                      |
|s_trajectories_c_mod.fpp |s_trajectory_c 初期化/コピー/フレーム入手/メモリ確保/解放など|

### 共通プログラム関連

|ファイル名               |内容                                           |
|:------------------------|:----------------------------------------------|
|c2py_util.py             |C データ -> numpy NDArray 変換                 |
|py2c_util.py             |numpy NDArray -> C データ変換                  |
|conv_f_c_util.fpp        |C データ <-> FORTRAN データ変換                |
|genesis_exe.py           |Python から FORTRAN 関数呼び出す               |
|libgenesis.py            |FORTRAN サブルーチンの C 型定義                |
|mbar_rename.py           |mbar_analysis 関連の module 名置換ファイル生成ツール|
|ctrl_c_mod.fpp           |iso_c_binding版コントロールデータ              |
|ctrl_files.py            |一時的コントロールファイル出力                 |

### crd_convert 関連

|ファイル名               |内容                                           |
|:------------------------|:----------------------------------------------|
|crd_convert.py           |crd_convert 回帰テスト実行プログラム           |
|crd_convert_c_mod.fpp    |GENESIS の cc_main.fpp と cc_setup.fpp の機能  |
|crd_convert_convert.fpp  |GENESIS の cc_convert.fpp の機能               |

### trj_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|trj_analysis.py           |trj_analysis 回帰テスト実行プログラム(Distance/Angle/Dihedral)|
|trj_analysis_c_mod.fpp    |GENESIS の ta_main.fpp と ta_setup.fpp の機能|
|trj_analysis_analysis.fpp |GENESIS の ta_analyze.fpp の機能             |

### wham_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|wham_analysis.py          |wham_analysis 回帰テスト実行プログラム       |
|wa_analysis_c_mod.fpp     |GENESIS の wa_main.fpp と wa_setup.fpp の機能|
|wa_analysis_analysis.fpp  |GENESIS の wa_analyze.fpp の機能             |

### mbar_analysis 関連

|ファイル名                |内容                                                        |
|:-------------------------------|:-----------------------------------------------------|
|mbar_analysis.py                |mbar_analysis 回帰テスト実行プログラ                  |
|mbar_analysis_umbrella_1d.py    |mbar_analysis 回帰テスト実行プログラム(Umbrella 1D)   |
|mbar_analysis_umbrella_block.py |mbar_analysis 回帰テスト実行プログラム(Umbrella Block)|
|mbar_analysis_c_mod.fpp         |GENESIS の ma_main.fpp と ma_setup.fpp の機能         |
|mbar_analysis_analysis.fpp      |GENESIS の ma_analyze.fpp の機能                      |

* mbar_rename.py 実行後に genesis/src/analysis/interface/mbar_analysis に mbar_analysis 関連の module 名変更ファイルが存在

### avecrd_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|avecrd_analysis.py        |avecrd_analysis 回帰テスト実行プログラム     |
|aa_analysis_c_mod.fpp     |GENESIS の aa_main.fpp と aa_setup.fpp の機能|
|aa_analysis_analysis.fpp  |GENESIS の aa_analyze.fpp の機能             |

### kmeans_clustering 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|kmeans_clustering.py      |kmeans_clustering 回帰テスト実行プログラム   |
|kc_analysis_c_mod.fpp     |GENESIS の kc_main.fpp と kc_setup.fpp の機能|
|kc_analysis_analysis.fpp  |GENESIS の kc_analyze.fpp の機能             |

### hb_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|hb_analysis_count_atom.py |hb_analysis 回帰テスト実行プログラム         |
|hb_analysis_count_snap.py |hb_analysis 回帰テスト実行プログラム         |
|hb_analysis_c_mod.fpp     |GENESIS の hb_main.fpp と hb_setup.fpp の機能|
|hb_analysis_analysis.fpp  |GENESIS の hb_analyze.fpp の機能             |

### rmsd_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|rmsd_analysis.py          |rmsd_analysis 回帰テスト実行プログラム       |
|ra_analysis_c_mod.fpp     |GENESIS の ra_main.fpp と ra_setup.fpp の機能|
|ra_analysis_analysis.fpp  |GENESIS の ra_analyze.fpp の機能             |

### drms_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|drms_analysis.py          |drms_analysis 回帰テスト実行プログラム       |
|dr_analysis_c_mod.fpp     |GENESIS の dr_main.fpp と dr_setup.fpp の機能|
|dr_analysis_analysis.fpp  |GENESIS の dr_analyze.fpp の機能             |

### rg_analysis 関連

|ファイル名                |内容                                         |
|:-------------------------|:--------------------------------------------|
|rg_analysis.py            |rg_analysis 回帰テスト実行プログラム         |
|rg_analysis_c_mod.fpp     |GENESIS の rg_main.fpp と rg_setup.fpp の機能|
|rg_analysis_analysis.fpp  |GENESIS の rg_analyze.fpp の機能             |

### msd_analysis 関連

|ファイル名                  |内容                                         |
|:---------------------------|:--------------------------------------------|
|msd_analysis.py             |msd_analysis 回帰テスト実行プログラム         |
|ma_analysis_c_mod.fpp       |GENESIS の ma_main.fpp と ma_setup.fpp の機能|
|ma_analysis_analysis.fpp    |GENESIS の ma_analyze.fpp の機能             |

### diffusion_analysis 関連

|ファイル名                        |内容                                         |
|:---------------------------------|:--------------------------------------------|
|diffusion_analysis.py             |diffusion_analysis 回帰テスト実行プログラム  |
|diffusion_analysis_main_c_mod.fpp |GENESIS の da_main.fpp と da_setup.fpp の機能|
|diffusion_analysis_analyze.fpp    |GENESIS の da_analyze.fpp の機能             |

### MDTraj, MDAnalysis 関連

|ファイル名                       |内容                                 |
|:---------------------------------|:-----------------------------------|
|test_mdanalysis.py                |MdAnalysis 回帰テスト実行プログラム |
|test_mdtraj.py                    |MdTraj 回帰テスト実行プログラム     |