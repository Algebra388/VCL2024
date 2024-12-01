## Task1

根据 `phong` 光照和 `blinn phong` 光照的公式即可完成渲染。
`phong` 为 $k_d(I_a + I_d \max(0 , n\dot l)) + k_s I_d \max(0, dir\dot r)^p$ 其中 $r$ 为反射光线， $dir$ 为视角方向。
`blinn phong` 为 $k_d(I_a + I_d \max(0 , n\dot l)) + k_s I_d \max(0, n\dot h)^p$ 其中 $h$ 为半程向量。

blinn phong 效果

![](~/home/algebra/zlt_codes/Grade2.1/vcl2024/vcx2024/build/linux/x86_64/release/1blinnphong.png)

phong 效果

![](~/home/algebra/zlt_codes/Grade2.1/vcl2024/vcx2024/build/linux/x86_64/release/1phong.png)

1. 顶点着色器和片段着色器的关系是什么样的？顶点着色器中的输出变量是如何传递到片段着色器当中的？

顶点着色器将顶点属性处理之后传递给片段着色器，将顶点着色器的输出作为片段着色器的输入。顶点着色器处理顶点后将数据传递给片段着色器通过插值来渲染整个图形。

再顶点着色器声明一个与片段着色器 `in` 属性的变量名字相同的 `out` 变量，再链接的时候就能进行传递。

2. 代码中的 `if (diffuseFactor.a < .2) discard;` 这行语句，作用是什么？为什么不能用 `if (diffuseFactor.a == 0.) discard;` 代替？

作用是当这个片段的 `alpha` 值 $< 0.2$ 就不进行绘制，因为 `alpha` $<0.2$ 就几乎透明了，不然很透明的地方会被渲染得很暗。

## Task3

1.使用 glCullFace 函数，可以剔除几何体的正面或背面。
正面剔除（GL_FRONT）：渲染背面。
背面剔除（GL_BACK）：渲染正面。

2.线的宽度难以控制。

## Task4

1.有向光源:透视投影矩阵
  平行光源:正交矩阵
2.因为坐标已经到了 $[-1,1]$ ，通过 `pos=pos * 0.5 + 0.5` 可以算出对应深度。