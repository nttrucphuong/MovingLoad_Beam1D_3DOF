using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Constraint
    {
        private Node node;
        public bool[] IsFixed { get; set; }
        public Constraint(Node node, params bool[] isFixed)
        {
            this.node = node;
            IsFixed = isFixed;
        }
        public Node GetNode()
        {
            return node;
        }
    }
}
