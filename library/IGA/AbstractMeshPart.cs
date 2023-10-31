//chuyển qua NURBS
//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace DEMSoft.IGA
//{
//    public abstract class AbstractMeshPart
//    {
//        private int countField = -1;
//        private int countDimension = -1;
//        private int[] tArray;
//        private int[] tArrayGlobal;

//        public void SetNumberOfFields(int d)
//        {
//            countField = d;
//        }

//        public int GetNumberOfFields()
//        {
//            return countField;
//        }

//        public void SetDimension(int d)
//        {
//            countDimension = d;
//        }

//        public int CountDimension()
//        {
//            return countDimension;
//        }

//        public void SetTArray(int[] tArray)
//        {
//            this.tArray = tArray;
//        }

//        public virtual int[] GetTArray()
//        {
//            return tArray;
//        }

//        public void SetTArrayGlobal(int[] tArrayGlobal)
//        {
//            this.tArrayGlobal = tArrayGlobal;
//        }

//        public virtual int[] GetTArrayGlobal()
//        {
//            return tArrayGlobal;
//        }
//    }
//}
