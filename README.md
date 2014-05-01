# 音声分析合成方式 World + waf
------------------------------------

## これは何？

* worldをwafでビルド・インストールできるようにしたものです。
* Cからも呼べるようにextern "C" をつけました
* Goから呼びたかったので作りました
* 現在は、WORLD v0.1.4（2014/05/01 での最新は0.1.4）を元になっています

worldについては、[本家ホームページ](http://ml.cs.yamanashi.ac.jp/world/)を参照ください

## インストール

     ./waf configure && ./waf
     sudo ./waf install

## サンプル

### Go の場合

      go run world.go

### C++ の場合

      g++ world.cpp -lworld

### C の場合

      gcc world.c -lworld

## 注意

* 公式のパッケージではありません
* 変更が、作者の意図に沿っていない可能性があります

## ソースコードに加えた変更点

* 主な音声分析変換合成のAPIを提供するヘッダファイルに extern "C" を入れた
* 名前衝突を回避するため、同一名の関数のうちdeprecatedになるであろう関数名を変更した（e.g. Dio -> DioOld）

Goから使おうと思ったらこうなりました。おしまい
