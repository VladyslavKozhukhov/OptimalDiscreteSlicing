using g3;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ErrDic = System.Collections.Generic.Dictionary<System.Tuple<int, int>, System.Tuple<int, int>>;
using TupleDInt = System.Tuple<int, int>;
using TupleErr_Sum = System.Tuple<System.Collections.Generic.Dictionary<System.Tuple<int, int>, int>, System.Collections.Generic.Dictionary<System.Tuple<int, int, int, int>, int>>;

namespace OptimalDiscreteSlicing
{
    class Program
    {

        static bool test = false;
        static bool weightsOn = false;
        public static ErrDic optDiscreteSlicingAlgo(Dictionary<TupleDInt, int> errorDic, HashSet<int> height, int N)
        {
            ErrDic ETmy_t_err = new ErrDic();
            Dictionary<Tuple<int, int>, int> E = new Dictionary<TupleDInt, int>();
            Dictionary<Tuple<int, int>, int> phi = new Dictionary<Tuple<int, int>, int>();
            int tMax = height.Max();
            int tMin = height.Min();

            for (int y = 1 - tMax; y <= 0; y++)
            {
                TupleDInt tmp = new TupleDInt(0, y);
                E[tmp] = 0;
            }
            for (int y = 1; y <= N + tMax - 1; y++)
            {
                for (int m = 1; m <= (int)Math.Round((double)((N - 2) / tMin)) + 1; m++)
                {
                    TupleDInt tmp = new TupleDInt(m, y);
                    E[tmp] = int.MaxValue;//inf
                    foreach (int t in height)
                    {
                        TupleDInt tmpNew = new TupleDInt(m - 1, y - t);

                        TupleDInt tmp1 = new TupleDInt(y - t, y);
                        TupleDInt tmp2 = new TupleDInt(m, y);
                        var errFromDic = -1;
                        var errInit = -1;
                     
                        if (E.ContainsKey(tmpNew)) { errInit = E[tmpNew]; }
                        else errInit = int.MaxValue;

                        errFromDic = errorDic[tmp1];//134 136


                        if (errInit + errFromDic < E[tmp2] && errInit != int.MaxValue)
                        {
                            phi[tmp2] = t;
                            E[tmp2] = errInit + errFromDic;
                        }
                       
                    }
                }
            }            
            foreach (TupleDInt tpl in phi.Keys)
            {
                var value = new TupleDInt(phi[tpl], E[tpl]);
                ETmy_t_err[tpl] = value;
            }
            return ETmy_t_err;

        }
        public static Tuple<int, int> findStartPoint(Dictionary<Tuple<int, int>, Tuple<int, int>> ETmy_t_err, int N, int tMax, int tMin)
        {
            int y = N - 1;
            int minErr = -1;
            bool firstIter = true;
            int yMin = -1;
            int mMin = -1;
           
            foreach (TupleDInt kk in ETmy_t_err.Keys)
            {
                if (kk.Item2 >= y && kk.Item2 <= N + tMax - 1)
                {

                    if (firstIter)
                    {                     
                        minErr = ETmy_t_err[kk].Item2;
                        firstIter = false;
                        yMin = kk.Item2;
                        mMin = kk.Item1;                  
                    }
                    else
                    {
                        if (minErr > ETmy_t_err[kk].Item2)
                        {
                            minErr = ETmy_t_err[kk].Item2;
                            firstIter = false;
                            yMin = kk.Item2;
                            mMin = kk.Item1;
                        }
                    }
                }
            }
            Tuple<int, int> retVal = new Tuple<int, int>(mMin, yMin);

            return retVal;
        }


        public static int getOptMperConstY(Dictionary<Tuple<int, int>, Tuple<int, int>> ETmy_t_err, int y, int tMin, int N, Dictionary<int, double> weights)
        {
            int mOpt = -1;
            double minErr = -1;
            Boolean first = true;
            foreach (TupleDInt kk in ETmy_t_err.Keys)
            {
                if (!weightsOn)
                {
                    if (kk.Item2 == y)
                    {
                        if (first) { minErr = ETmy_t_err[kk].Item2; mOpt = kk.Item1; first = false; }                     
                        else if (minErr > ETmy_t_err[kk].Item2)
                        {
                            mOpt = kk.Item1;
                        }
                    }
                }
                else
                {

                    if (kk.Item2 == y)
                    {
                        if (first) { minErr = ETmy_t_err[kk].Item2*weights[ETmy_t_err[kk].Item1]; mOpt = kk.Item1; first = false; }                
                        else if (minErr > ETmy_t_err[kk].Item2 * weights[ETmy_t_err[kk].Item1])
                        {
                            mOpt = kk.Item1;
                        }
                    }
                }
            }
            return mOpt;
        }





        public static List<int> getOptSlice(Tuple<int, int> startPoint, Dictionary<Tuple<int, int>, Tuple<int, int>> ETmy_t_err, int tMin, int N, Dictionary<int, double> weights)
        {
            List<int> path = new List<int>();
            Tuple<int, int> tmp = new Tuple<int, int>(startPoint.Item1, startPoint.Item2);
            int yNew = startPoint.Item2;
            int mNew = -1;
            Boolean first = true;
            path.Add(startPoint.Item2);

            while (yNew > 1 || first)//Check if stop on 1 or 0!
            {

                yNew = yNew - ETmy_t_err[tmp].Item1; //z-t
                mNew = getOptMperConstY(ETmy_t_err, yNew, tMin, N,weights);
                path.Add(yNew);
                if (first) first = false; else first = false;
                tmp = new Tuple<int, int>(mNew, yNew);
            }
            return path;
        }


        public static List<int> getIntersections(int x, int z, Bitmap3 bmp)
        {
            List<int> interList = new List<int>();
            bool inShape = false;
            for (int y = 0; y < bmp.Dimensions.y; y++)
            {
                //we were out of shape and now going in shape
                if (bmp.Get(createVector(x, y, z)) == true && inShape == false)
                {
                    inShape = true;
                    interList.Add(y);
                }
                //we were in shape and now going out of shape
                if (bmp.Get(createVector(x, y, z)) == false && inShape == true)
                {
                    inShape = false;
                    interList.Add(y);
                }
                //obj reaches the pick
                if (bmp.Get(createVector(x, y, z)) == true && y == bmp.Dimensions.y - 1)
                {
                    interList.Add(y + 1);
                }
            }
            return interList;
        }

        public static Tuple<Dictionary<Tuple<int, int>, int>, Dictionary<Tuple<int, int, int, int>, int>> calcErrorAndSum(Bitmap3 bmp, int tmax, int tmin)
        {
            Dictionary<Tuple<int, int>, int> errorDic = new Dictionary<Tuple<int, int>, int>(); //key: zi,zj value: total error
            Dictionary<Tuple<int, int, int, int>, int> sumDic = new Dictionary<Tuple<int, int, int, int>, int>(); //key: zi,zj,x,z value: sum
            List<int> intersectionList;

            //zerro all errors
            for (int zi = 1 - tmax; zi <= bmp.Dimensions.y + tmax; zi++)
            {
                for (int zj = zi; zj <= Math.Min(zi + tmax, bmp.Dimensions.y + tmax); zj++)
                {
                    Tuple<int, int> key = new Tuple<int, int>(zi, zj);
                    errorDic[key] = 0;
                }
            }

            //run for each x and z (2D looking from top of object)
            for (int x = 0; x < bmp.Dimensions.x; x++)
            {
                for (int z = 0; z < bmp.Dimensions.z; z++)
                {
                    //calculate error and sum for specific x and y
                    intersectionList = getIntersections(x, z, bmp);
                    for (int k = 0; k < intersectionList.Count; k++)
                    {
                        //check for null pointer exception 
                        int zi;
                        if (k > 0)
                        {
                            zi = Math.Max(intersectionList[k] - tmax, intersectionList[k - 1]);
                        }
                        else
                        {
                            zi = intersectionList[k] - tmax;
                        }
                        for (; zi <= intersectionList[k]; zi++)
                        { //check ranges!
                            for (int zj = intersectionList[k]; zj <= zi + tmax; zj++)
                            {
                                int s = (int)Math.Pow(-1, k) * (zi - intersectionList[k]);
                                int l = k + 1;
                                while (l < intersectionList.Count && intersectionList[l] < zj)
                                {
                                    s += (int)Math.Pow(-1, l) * (intersectionList[l - 1] - intersectionList[l]);
                                    l++;
                                }
                                s += (int)Math.Pow(-1, l) * (intersectionList[l - 1] - zj);
                                Tuple<int, int, int, int> sumKey = new Tuple<int, int, int, int>(zi, zj, x, z);
                                Tuple<int, int> key = new Tuple<int, int>(zi, zj);
                                if (!sumDic.ContainsKey(sumKey)) //CHECK THIS OUT IT LOOKS BAD!!! (IN THE LOOK WE REPEAT CALCULATION FOR ZI,ZJ,X,Z!!!!
                                {
                                    sumDic.Add(sumKey, s);
                                    int error = zj - zi - s;
                                    if (errorDic.ContainsKey(key)) //check if need to add prev error
                                    {
                                        errorDic[key] = errorDic[key] + error;
                                    }
                                    else
                                    {
                                        errorDic[key] = error;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return new Tuple<Dictionary<Tuple<int, int>, int>, Dictionary<Tuple<int, int, int, int>, int>>(errorDic, sumDic);
        }

        public static Vector3i createVector(int x, int y, int z)
        {
            return new Vector3i(x, y, z);
        }

        public static Bitmap3 createVoxelizedRepresentation(String objPath)
        {
            DMesh3 mesh = StandardMeshReader.ReadMesh(objPath);

            int num_cells = 128;
            double cell_size = mesh.CachedBounds.MaxDim / num_cells;

            MeshSignedDistanceGrid sdf = new MeshSignedDistanceGrid(mesh, cell_size);
            sdf.Compute();

            //** voxels**//
            Bitmap3 bmp = new Bitmap3(sdf.Dimensions);
            Console.WriteLine(bmp.Dimensions.x + " " + bmp.Dimensions.y + " " + bmp.Dimensions.z);
            foreach (Vector3i idx in bmp.Indices())
            {
                float f = sdf[idx.x, idx.y, idx.z];
                bmp.Set(idx, (f < 0) ? true : false);
                //for bunny only removes bottom
                if (idx.y < 8)
                {
                    bmp.Set(idx, false);
                }

                if (test) //take only one line from top
                {
                    if (idx.z != 50 || idx.x != 60)
                    {
                        bmp.Set(idx, false);
                    }
                    else
                    {
                        bmp.Set(idx, true);
                        Console.WriteLine(bmp.Get(idx));
                    }
                }

            }
            return bmp;
        }

        public static void printVoxelizedRepresentation(Bitmap3 bmp, String outputPath)
        {
            VoxelSurfaceGenerator voxGen = new VoxelSurfaceGenerator();
            voxGen.Voxels = bmp;      
            voxGen.Generate();
            DMesh3 voxMesh = voxGen.Meshes[0];
            Util.WriteDebugMesh(voxMesh, outputPath);
        }

        public static bool isCoulumnInObj(Bitmap3 bmp, int zi, int zj, int x, int z)
        {
            bool visitedInObject = false;
            for (int y = zi; y < zj; y++)
            {
                if (y < 0)
                {
                    continue;
                }
                if (y >= bmp.Dimensions.y)
                {
                    if(visitedInObject){
                         return true;
                    }
                    else{
                        return false;
                    }
                }
                if (!bmp.Get(createVector(x,y, z)))
                {
                    return false;
                }
                else
                {
                    visitedInObject = true;
                }
            }
            if (visitedInObject)
            {
                return true;
            }
            return false;
        }

        public static Bitmap3 createNewObjectForPriniting(List<int> path, Dictionary<Tuple<int, int, int, int>, int> sumDic, Vector3i newObjDim, Bitmap3 oldObj)
        {
            Bitmap3 voxPrintResult = new Bitmap3(newObjDim);
            foreach (Vector3i idx in voxPrintResult.Indices()) //initialize object
            {
                voxPrintResult.Set(idx, false);
            }
            for (int i = 0; i < path.Count - 1; i++)
            { //-1 because we want to take every time the range Zj to Zi where Zj > Zi
                int zj = path[i]-path.Last();
                int zi = path[i + 1] - path.Last();
                for (int x = 0; x < voxPrintResult.Dimensions.x; x++)
                {
                    for (int z = 0; z < voxPrintResult.Dimensions.z; z++)
                    {
                        Tuple<int, int, int, int> key = new Tuple<int, int, int, int>(zi, zj, x, z);
                        if ((sumDic.ContainsKey(key) && sumDic[key] >= 0) || isCoulumnInObj(oldObj, zi, zj, x, z))
                        {
                            for (int y = zi; y < zj; y++)
                            {
                                voxPrintResult.Set(createVector(x, y, z), true);
                            }
                        }
                    }
                }
            }
            return voxPrintResult;
        }
        

        static void Main(string[] args)
        {
          
            Console.WriteLine("Path: ");
            string pathToFile = Console.ReadLine();
            Console.WriteLine("Weight Mode?(Y/N): ");
            string answer = Console.ReadLine();
            List<int> input= new List<int>();
            Dictionary<int, double> weights = new Dictionary<int, double>();
            if (answer.ToUpper() == "Y")
            {
                int size = -1;
                int enter = 1;  
                while (enter>0) { 
                    Console.WriteLine("Please enter slice size: ");
                    size = Convert.ToInt32(Console.ReadLine());
                    if(size <= 0)
                    {
                        enter = -1;
                    }
                    if (enter > 0)
                    {
                        Console.WriteLine("weight : ");
                        if (!weights.ContainsKey(size) && enter > 0)
                        {
                            double wg = -1;
                            wg = Convert.ToDouble(Console.ReadLine());
                            if (wg > 0.0) {
                                weights.Add(size, wg);
                                input.Add(size);
                            }
                            else
                            {
                                Console.WriteLine("set default weight 1");
                                weights.Add(size, 1);
                                input.Add(size);
                            }
                        }
                        
                    }                    

                }
                weightsOn = true;
            }
            else
            {
                Console.WriteLine("Insert heights: ");
                string s = Console.ReadLine();
                input = s.Split(' ').Select(t => Convert.ToInt32(t)).ToList<int>();
                weightsOn =false;
            }
         

            HashSet<int> legitSliceHights = new HashSet<int>(input);
            //"C:\\Users\\VladKo\\Downloads\\bunny.obj"
            Bitmap3 bmp = createVoxelizedRepresentation(pathToFile);
            printVoxelizedRepresentation(bmp, "C:\\Users\\VladKo\\Downloads\\inputVox.obj");
            if (test)
            {
                getIntersections(60, 50, bmp).ForEach(Console.WriteLine);
            }
            Tuple<Dictionary<Tuple<int, int>, int>, Dictionary<Tuple<int, int, int, int>, int>> errorAndSum = calcErrorAndSum(bmp, legitSliceHights.Max(), legitSliceHights.Min());
            Dictionary<Tuple<int, int>, Tuple<int, int>> algResults = optDiscreteSlicingAlgo(errorAndSum.Item1, legitSliceHights, bmp.Dimensions.y);
            Tuple<int, int> startPoint = findStartPoint(algResults, bmp.Dimensions.y, legitSliceHights.Max(), legitSliceHights.Min());
            List<int> path = getOptSlice(startPoint, algResults, legitSliceHights.Min(), bmp.Dimensions.y,weights); //from top to bottom
            Vector3i newObjDim = createVector(bmp.Dimensions.x, path.First() - path.Last(), bmp.Dimensions.z);
            Bitmap3 outputObj = createNewObjectForPriniting(path, errorAndSum.Item2, newObjDim, bmp);
            printVoxelizedRepresentation(outputObj, "C:\\Users\\VladKo\\Downloads\\outputVox.obj");

        }
    }
}