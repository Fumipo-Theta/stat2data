# 一様分布から値を取り出す関数. 値域が2つの領域に分かれている場合にも対応する.
def randomFrom(range1,range2 = []):
    if (len(range2) > 2):
        f = (range1[1]-range1[0] + range2[1]-range2[0])/(range1[1]-range1[0])
        r = rand()
        if (r > f):
            return (range2[1]-range2[0])*(r-f)/(1-f) + range2[0]
        else:
            return (range1[1]-range1[0])*r/f + range1[0]
        
    else:
        return (range1[1]-range1[0]) * rand() + range1[0]

# 半径1の球面上の座標を取得する. (X,Y,Z)に対応する.
def onSuperSphere(mu,sigma,minmax):
    R = np.sqrt(3)*sigma
    u = 2 * (minmax[1] - mu) # 変数の上限による
    l = 2 * (mu - minmax[0]) # 変数の下限による
    
    # 四面体による断面の小円の半径の二乗
    dl = R * R - l * l /3
    du = R * R - u * u /3
    
    if ((np.sqrt(3) * l < R) | (np.sqrt(3) * u < R) ):
        # 断面の小円の半径が球の半径と一致するとき球面はなくなる
        xyz = noSphere()
    
    elif ((dl >= 0) & (du >=0)):
        # 下限と上限の両方に対応する小円が球面に現れるとき
        xyz = cutByOctahedron(l,u,R)
        
    elif (dl >= 0):
        # 下限の小円のみ現れるとき
        xyz = cutByLowerTetrahedron(l,u,R)
    elif (du >= 0):
        # 上限の小円のみ現れるとき
        xyz = cutByUpperTetrahedron(l,u,R)
    else:
        # 四面体による断面が現れないとき
        xyz = completeSphere(l,u,R)
    
    return xyz

# 球面が2つの正四面体で切られる場合
def cutByOctahedron(m,l,R):
    pi = np.pi
    
    # 切り口の小円の半径
    rm = np.sqrt(R*R - m*m/3)
    rl = np.sqrt(R*R - l*l/3)
    
    # 小円のZ座標の上限下限
    zm_top = (m + np.sqrt(6) * rm)/3
    zm_bottom = (m - np.sqrt(6) * rm)/3
    
    zl_top = (l + np.sqrt(6) * rl)/3
    zl_bottom = (l - np.sqrt(6) * rl)/3
    
    # 下限と上限の小円の半径の内どちらが大きいかで条件が異なる
    if (rm >= rl):
        # 最小値の小円のほうが大きいとき
        
        if (np.sqrt(2) * R <= m ):
            # 小円が半球内に収まるとき
            z_max = R
            z_min = 0
            z = randomFrom([z_min,z_max])
            r = np.sqrt(R*R - z*z)
            
            bm = np.sqrt(2)/2 * (m-z)
            bl = np.sqrt(2)/2 * (l-z)
            
            if (zm_top <= z):
                theta = randomFrom(
                    [
                        pi/4,
                        pi*3/4
                    ])
                
            elif (zl_top <= z ):
                theta = randomFrom(
                    [
                        pi/4 + np.arccos(bm/r),
                        pi*3/4
                    ])
                
            elif (zl_bottom <= z):
                theta = randomFrom(
                    [
                        np.arccos(bm/r)+pi/4,
                        pi*3/4 - np.arccos(bl/r)
                    ])
                
            elif (zm_bottom <= z):
                theta = randomFrom(
                [
                    pi/4 + np.arccos(bm/r),
                    pi*3/4
                ])
                
            else:
                theta = randomFrom([pi/4,pi*3/4])
            
        else: 
            # 小円が赤道を切るとき
            
            if(R <= m):
                # 小円mどうしは交差しないとき
                z_max = R
                z_min = 0
            else:
                # 同じタイプの小円が交差するとき
                z_max = m
                z_min = np.sqrt((R*R-m*m)/2) # 小円の交差点
            
            # 最大値と最小値に対応する2種類の小円が交差するかどうか
            d_intersect_double_circles = 8 * R*R + 2 * m * l - 3 * (m*m + l*l)
            
            if (d_intersect_double_circles > 0):
                # 2タイプの小円が交差するとき
                z_intersect_top =    (m+l + np.sqrt(d_intersect_double_circles))/4
                z_intersect_bottom = (m+l - np.sqrt(d_intersect_double_circles))/4

                # もし小円mの交点より上に小円mとlの交点がある場合は解なし
                if (z_max < z_intersect_top):
                    return {
                        "x":np.nan,
                        "y":np.nan,
                        "z":np.nan
                    }
                
                
                z = randomFrom(
                    [z_min,z_intersect_bottom],
                    [z_intersect_top,z_max]
                )
                r = np.sqrt(R*R - z*z)
            
                bm = np.sqrt(2)/2 * (m-z)
                bm_inv = np.sqrt(2)/2 * (m+z)
                bl = np.sqrt(2)/2 * (l-z)
                
                
                
                # 2タイプの小円が交差する
                if (z >= zm_top):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4
                        ]
                        
                    )
                elif (z >= zl_top):
                    theta = randomFrom(
                        [
                            pi/4 + np.arccos(bm/r),
                            pi*3/4
                        ]
                    )
                
                elif (z >= z_intersect_top):
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )

                elif (z >= (l-m)/2):
                    # 小円lと, 南半球の小円mの交点
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )

                else:
                    # 南半球の小円mとの間
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bm_inv/r)
                        ]
                    )
                
            else:
                # 2タイプの小円は交差しない
                z = randomFrom(
                    [z_min,z_max],
                )
                r = np.sqrt(R*R - z*z)
            
                bm = np.sqrt(2)/2 * (m-z)
                bm_inv = np.sqrt(2)/2 * (m+z)
                bl = np.sqrt(2)/2 * (l-z)
                
                if (zm_top < z):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4
                    ])
                elif (zl_top < z):
                    theta = randomFrom(
                        [
                            pi/4 + np.arccos(bm/r),
                            pi*3/4
                    ])
                elif (zl_bottom < z):
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )
                elif (-zm_bottom < z):
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4
                        ]
                    )
                else:
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bm_inv/r)
                        ]
                    )
            
    
    else:
        # 小円lのほうが大きいとき
        
        if (np.sqrt(2) * R <= l ):
            # 小円が半球内に収まるとき
            z_max = R
            z_min = 0
            z = randomFrom([z_min,z_max])
            r = np.sqrt(R*R - z*z)
            
            bm = np.sqrt(2)/2 * (m-z)
            bl = np.sqrt(2)/2 * (l-z)
            
 
            
            if (zl_top <= z):
                theta = randomFrom(
                    [
                        pi/4,
                        pi*3/4
                    ])
                
            elif (zm_top <= z ):
                theta = randomFrom(
                    [
                        pi/4 ,
                        pi*3/4 - np.arccos(bl/r)
                    ])
                
            elif (zm_bottom <= z):
                theta = randomFrom(
                    [
                        np.arccos(bm/r)+pi/4,
                        pi*3/4 - np.arccos(bl/r)
                    ])
                
            elif (zl_bottom <= z):
                theta = randomFrom(
                [
                    pi/4 ,
                    pi*3/4 - np.arccos(bm/r)
                ])
                
            else:
                theta = randomFrom([pi/4,pi*3/4])
            
        else: 
            # 小円が赤道を切るとき
            
            if(R <= l):
                # 小円mどうしは交差しないとき
                z_max = R
                z_min = 0
            else:
                # 同じタイプの小円が交差するとき
                z_max = l
                z_min = np.sqrt((R*R-l*l)/2) # 小円の交差点
            
            
            
            # 最大値と最小値に対応する2種類の小円が交差するかどうか
            d_intersect_double_circles = 8 * R*R + 2 * m * l - 3 * (m*m + l*l)

            if (d_intersect_double_circles > 0):
                # 2タイプの小円が交差するとき
                z_intersect_top = (m+l + np.sqrt(d_intersect_double_circles))/4
                z_intersect_bottom = (m+l - np.sqrt(d_intersect_double_circles))/4

                # もし小円mの交点より上に小円mとlの交点がある場合は解なし
                if (z_max < z_intersect_top):
                    return {
                        "x":np.nan,
                        "y":np.nan,
                        "z":np.nan
                    }
                
                
                z = randomFrom(
                    [z_min,z_intersect_bottom],
                    [z_intersect_top,z_max]
                )
                r = np.sqrt(R*R - z*z)
            
                bm = np.sqrt(2)/2 * (m-z)
                bl_inv = np.sqrt(2)/2 * (l+z)
                bl = np.sqrt(2)/2 * (l-z)
                
                
                
                # 2タイプの小円が交差する
                if (z >= zl_top):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4
                    ])
                elif (z >= zm_top):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4- np.arccos(bl/r)
                    ])
                
                elif (z >= z_intersect_top):
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )

                elif (z >= (m-l)/2):
                    # 小円lと, 南半球の小円mの交点
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )

                else:
                    # 南半球の小円mとの間
                    theta = randomFrom(
                        [
                            np.arccos(bl_inv/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )
                
            else:
                # 2タイプの小円は交差しない
                z = randomFrom(
                    [z_min,z_max],
                )
                r = np.sqrt(R*R - z*z)
            
                bm = np.sqrt(2)/2 * (m-z)
                bl_inv = np.sqrt(2)/2 * (l+z)
                bl = np.sqrt(2)/2 * (l-z)
                
                if (zl_top < z):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4
                    ])
                elif (zm_top < z):
                    theta = randomFrom(
                        [
                            pi/4 ,
                            pi*3/4-np.arccos(bl/r)
                    ])
                elif (zm_bottom < z):
                    theta = randomFrom(
                        [
                            np.arccos(bm/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )
                elif (-zl_bottom < z):
                    theta = randomFrom(
                        [
                            pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )
                else:
                    theta = randomFrom(
                        [
                            np.arccos(bl_inv/r)+pi/4,
                            pi*3/4 - np.arccos(bl/r)
                        ]
                    )
    
    rotate = rand()
    if (rotate <= 0.25):
        theta = pi/2 - theta

    elif (rotate <= 0.5):
        theta = pi*3/2 - theta
    elif (rotate <= 0.75):
        theta = pi + theta
    
    if (rand() < 0.5):
        return {
            "x" : r * np.cos(theta),
            "y" : r * np.sin(theta),
            "z" : z
        }
    else:
        return {
            "x" : -r * np.cos(theta-pi/2),
            "y" : -r * np.sin(theta-pi/2),
            "z" : -z
        }

    
def cutByLowerTetrahedron(m,l,R):
    pi = np.pi
    
    # 切り口の小円の半径
    rm = np.sqrt(R*R - m*m/3)
    
    # 小円のZ座標の上限下限
    zm_top = (m + np.sqrt(6) * rm)/3
    zm_bottom = (m - np.sqrt(6) * rm)/3

    
    if (np.sqrt(2) * R <= m ):
        # 小円が半球内に収まるとき
        z_max = R
        z_min = 0
        z = randomFrom([z_min,z_max])
        r = np.sqrt(R*R - z*z)
            
        bm = np.sqrt(2)/2 * (m-z)

        if (zm_top <= z):
            theta = randomFrom(
                [
                    pi/4,
                    pi*3/4
                ])
                
            
        elif (zm_bottom <= z):
            theta = randomFrom([
                 pi/4 + np.arccos(bm/r),
                 pi*3/4
            ])
                
        else:
            theta = randomFrom([pi/4,pi*3/4])
            
    else: 
        # 小円が赤道を切るとき
            
        if(R <= m):
            # 小円mどうしは交差しないとき
            z_max = R
            z_min = 0
        else:
            # 同じタイプの小円が交差するとき
            z_max = m
            z_min = np.sqrt((R*R-m*m)/2) # 小円の交差点
            
        # 2タイプの小円は交差しない
        z = randomFrom(
            [z_min,z_max],
        )
        r = np.sqrt(R*R - z*z)
            
        bm = np.sqrt(2)/2 * (m-z)
        bm_inv = np.sqrt(2)/2 * (m+z)

                
        if (zm_top < z):
            theta = randomFrom(
                [
                    pi/4,
                    pi*3/4
                ])

        elif (-zm_bottom < z):
            theta = randomFrom(
                [
                    pi/4 + np.arccos(bm/r),
                    pi*3/4
                ]
            )
        else:
            theta = randomFrom(
                [
                    np.arccos(bm/r)+pi/4,
                    pi*3/4 - np.arccos(bm_inv/r)
                ]
            )
    
    rotate = rand()
    if (rotate <= 0.25):
        theta = pi/2 - theta

    elif (rotate <= 0.5):
        theta = pi*3/2 - theta
    elif (rotate <= 0.75):
        theta = pi + theta
    
    if (rand() < 0.5):
        return {
            "x" : r * np.cos(theta),
            "y" : r * np.sin(theta),
            "z" : z
        }
    else:
        return {
            "x" : -r * np.cos(theta-pi/2),
            "y" : -r * np.sin(theta-pi/2),
            "z" : -z
        }
    
def cutByUpperTetrahedron(m,l,R):
    pi = np.pi
    
    # 切り口の小円の半径
    rl = np.sqrt(R*R - l*l/3)
    
    # 小円のZ座標の上限下限
    zl_top = (l + np.sqrt(6) * rl)/3
    zl_bottom = (l - np.sqrt(6) * rl)/3

    
    if (np.sqrt(2) * R <= l ):
        # 小円が半球内に収まるとき
        z_max = R
        z_min = 0
        z = randomFrom([z_min,z_max])
        r = np.sqrt(R*R - z*z)
            
        bl = np.sqrt(2)/2 * (l-z)

            
 
            
        if (zl_top <= z):
            theta = randomFrom(
                [
                    pi/4,
                    pi*3/4
                ])
                
            
        elif (zl_bottom <= z):
            theta = randomFrom([
                 pi/4 ,
                 pi*3/4 - np.arccos(bl/r)
            ])
                
        else:
            theta = randomFrom([pi/4,pi*3/4])
            
    else: 
        # 小円が赤道を切るとき
            
        if(R <= l):
            # 小円mどうしは交差しないとき
            z_max = R
            z_min = 0
        else:
            # 同じタイプの小円が交差するとき
            z_max = l
            z_min = np.sqrt((R*R-l*l)/2) # 小円の交差点
            
        # 2タイプの小円は交差しない
        z = randomFrom(
            [z_min,z_max],
        )
        r = np.sqrt(R*R - z*z)
            
        bl = np.sqrt(2)/2 * (l-z)
        bl_inv = np.sqrt(2)/2 * (l+z)

                
        if (zl_top < z):
            theta = randomFrom(
                [
                    pi/4,
                    pi*3/4
                ])

        elif (-zl_bottom < z):
            theta = randomFrom(
                [
                    pi/4,
                    pi*3/4 - np.arccos(bl/r),
                ]
            )
        else:
            theta = randomFrom(
                [
                    np.arccos(bl_inv/r)+pi/4,
                    pi*3/4 - np.arccos(bl/r)
                ]
            )
    
    rotate = rand()
    if (rotate <= 0.25):
        theta = pi/2 - theta

    elif (rotate <= 0.5):
        theta = pi*3/2 - theta
    elif (rotate <= 0.75):
        theta = pi + theta
    
    if (rand() < 0.5):
        return {
            "x" : r * np.cos(theta),
            "y" : r * np.sin(theta),
            "z" : z
        }
    else:
        return {
            "x" : -r * np.cos(theta-pi/2),
            "y" : -r * np.sin(theta-pi/2),
            "z" : -z
        }

    
def completeSphere(m,l,R):
    pi = np.pi
    z_max = R
    z_min = 0
    z = randomFrom([z_min,z_max])
    r = np.sqrt(R*R - z*z)
    
    theta = randomFrom(
                [
                    pi/4,
                    pi*3/4
                ])
    
    rotate = rand()
    if (rotate <= 0.25):
        theta = pi/2 - theta

    elif (rotate <= 0.5):
        theta = pi*3/2 - theta
    elif (rotate <= 0.75):
        theta = pi + theta
    
    if (rand() < 0.5):
        return {
            "x" : r * np.cos(theta),
            "y" : r * np.sin(theta),
            "z" : z
        }
    else:
        return {
            "x" : -r * np.cos(theta-pi/2),
            "y" : -r * np.sin(theta-pi/2),
            "z" : -z
        }

def noSphere():
    return {
        "x" : np.nan,
        "y" : np.nan,
        "z" : np.nan
    }