using DEMSoft.Drawing;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace OOFEM
{
    internal static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Geometry
            double L = 8; //beam's lenght
            double thickness = 1;
            double q = -1;
            double F = -10;
            double k = 1e6;

            //Material
            double E = 100;
            double nu = 0.25;
            double rho = 1;


            Console.WriteLine(theta[1]);
            Console.ReadKey();


            List<Node> nodes = new List<Node>();
            List<Node> nodes1 = new List<Node>();

            for (int i = 1; i <= numr + 1; i++)
            {
                double edge = (i - 1) * (R - r) / numr + r;
                for (int j = 1; j <= numc + 1; j++)
                {
                    nn += 1;
                    nodes.Add(new Node(Math.Cos(Math.PI / (2 * numc) * (j - 1)) * edge, Math.Sin(Math.PI / (2 * numc) * (j - 1)) * edge, 0));
                    nodes1.Add(new Node(Math.Cos(Math.PI / (2 * numc) * (j - 1)) * edge, Math.Sin(Math.PI / (2 * numc) * (j - 1)) * edge, 0));
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int numr2 = 5;
            int numc2 = 12;
            double R2 = 5;
            double r2 = 2;
            int nn2 = 0;

            List<Node> nodes2 = new List<Node>();

            for (int i = 1; i <= numr2 + 1; i++)
            {
                double edge = (i - 1) * (R2 - r2) / numr2 + r2;
                for (int j = 1; j <= numc2 + 1; j++)
                {
                    nn2 += 1;
                    nodes.Add(new Node(Math.Cos(Math.PI / (2 * numc2) * (j - 1)) * edge, Math.Sin(Math.PI / (2 * numc2) * (j - 1)) * edge, 0));
                    nodes2.Add(new Node(Math.Cos(Math.PI / (2 * numc2) * (j - 1)) * edge, Math.Sin(Math.PI / (2 * numc2) * (j - 1)) * edge, 0));
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            int[,] connectionMatrix = new int[numr * numc, 4];
            int c = 0;
            for (int i = 0; i < numr; i++)
            {
                for (int j = 0; j < numc; j++)
                {
                    connectionMatrix[c, 0] = (numc + 1) * i + j;
                    connectionMatrix[c, 1] = (numc + 1) * i + j + 1;
                    connectionMatrix[c, 2] = (numc + 1) * (i + 1) + j + 1;
                    connectionMatrix[c, 3] = (numc + 1) * (i + 1) + j;
                    c++;
                }
            }

            int[,] connectionMatrix2 = new int[numr2 * numc2, 4];
            int c2 = 0;

            for (int i = 0; i < numr2; i++)
            {
                for (int j = 0; j < numc2; j++)
                {
                    connectionMatrix2[c2, 0] = (numc2 + 1) * i + j;
                    connectionMatrix2[c2, 1] = (numc2 + 1) * i + j + 1;
                    connectionMatrix2[c2, 2] = (numc2 + 1) * (i + 1) + j + 1;
                    connectionMatrix2[c2, 3] = (numc2 + 1) * (i + 1) + j;
                    c2++;
                }
            }

            List<Element> elems = new List<Element>();
            List<Element> elems1 = new List<Element>();
            List<Element> elems2 = new List<Element>();

            for (int i = 0; i < numr; i++)
            {
                for (int j = 0; j < numc; j++)
                {
                    Plane2DElement element1 = new Plane2DElement(nodes1[connectionMatrix[numc * i + j, 0]],
                                                             nodes1[connectionMatrix[numc * i + j, 1]],
                                                             nodes1[connectionMatrix[numc * i + j, 2]],
                                                             nodes1[connectionMatrix[numc * i + j, 3]],
                                                             E,
                                                             nu,
                                                             thickness);
                    elems.Add(element1);
                    elems1.Add(element1);
                }
            }
            for (int i = 0; i < numr2; i++)
            {
                for (int j = 0; j < numc2; j++)
                {
                    Plane2DElement element1 = new Plane2DElement(nodes2[connectionMatrix2[numc2 * i + j, 0]],
                                                             nodes2[connectionMatrix2[numc2 * i + j, 1]],
                                                             nodes2[connectionMatrix2[numc2 * i + j, 2]],
                                                             nodes2[connectionMatrix2[numc2 * i + j, 3]],
                                                             E,
                                                             nu,
                                                             thickness);
                    elems.Add(element1);
                    elems2.Add(element1);
                }
            }


            List<Subdomain> subdomains = new List<Subdomain>();
            Subdomain subdomain1 = new Subdomain(nodes1, elems1);
            Subdomain subdomain2 = new Subdomain(nodes2, elems2);
            subdomains.Add(subdomain1);
            subdomains.Add(subdomain2);

            // interface
            List<Node> node_interface_subdomain1 = new List<Node>();
            List<Node> node_interface_subdomain2 = new List<Node>();

            for (int i = 0; i < numc + 1; i++)
            {
                node_interface_subdomain1.Add(nodes1[(numc + 1) * numr + i]);
            }

            for (int i = 0; i < numc2 + 1; i++)
            {
                node_interface_subdomain2.Add(nodes2[i]);
            }

            List<double[]> listContact = new List<double[]>();
            int numContact = 10;
            double lengthContact = 2;
            double dy = lengthContact / numContact;
            for (int i = 0; i < numContact; i++)
            {
                double xi = i * dy;
                double xin = (i + 1) * dy;
                listContact.Add(new double[] { xi, xin });
            }

            List<Node> listContact1 = new List<Node>();
            int numContact1 = 10;
            double R_contact = 2;
            for (int i = 1; i < numContact + 1; i++)
            {
                listContact1.Add(new Node(Math.Cos(Math.PI / (2 * numContact1) * (i - 1)) * R_contact, Math.Sin(Math.PI / (2 * numContact1) * (i - 1)) * R_contact, 0));
            }

            int numberofField = 2;
            StaticModel model = new StaticModel(numberofField);

            for (int i = 0; i < subdomains.Count; i++)
            {
                model.InsertSubdomain(subdomains[i]);
            }

            for (int i = 0; i < nodes.Count; i++)
            {
                model.InsertNode(nodes[i]);
            }

            for (int i = 0; i < elems.Count; i++)
            {
                model.InsertElement(elems[i]);
            }

            MortarContactElememt mortarContact = new MortarContactElememt(model.GetSubdomain(0), model.GetSubdomain(1), node_interface_subdomain1, node_interface_subdomain2, listContact);
            model.AddMortarContact(mortarContact);

            for (int j = 0; j < nodes.Count; j++)
            {
                if (nodes[j].GetPosition(1) == 0)
                {
                    Constraint cs = new Constraint(model.GetNode(j), false, true);
                    model.AddConstraint(cs);
                }
                else if (nodes[j].GetPosition(0) < 1e-6)
                {
                    Constraint cs = new Constraint(model.GetNode(j), true, false);
                    model.AddConstraint(cs);
                }
            }



            double q = thickness * p * r * Math.PI / (2 * numc);
            for (int j = 0; j < numc + 1; j++)
            {
                double fx = q * Math.Cos(Math.PI * j / (2 * numc));
                double fy = q * Math.Sin(Math.PI * j / (2 * numc));
                if (j == 0 || j == numc)
                {
                    fx = fx / 2.0;
                    fy = fy / 2.0;
                }
                Force f = new Force(nodes1[j], fx, fy);
                model.AddForce(f);
            }

            ViewerForm viewer = new ViewerForm(true);
            //model.DrawNodes(viewer);
            //model.DrawElements(viewer);
            //model.DrawConstrain(viewer);
            //model.DrawForce(viewer);


            //model.PreProcessing();
            //model.Solve();
            //model.PostProcessing();

            double scale = 20;
            //model.DrawDefomation(viewer, scale);
            model.DrawResult(viewer, scale, Result.USUM);


            //viewer.UpdateCamera();
            //viewer.Run();

        }
    }
}
