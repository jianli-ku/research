using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;
using TMPro;

//GlobalManager is for storing variables that will be set in the main menu and accessed by ARManager for calculation
public class GlobalManager : Singleton<GlobalManager>
{
    #region Private Variables

    //NOTE: These may be hard-set in the editor. You may have to manually go change their initial values there. 
    private static double lengthOption1 = 1.772f; //Note: when changing these values, be sure to adjust the displayed values in the dropdown menu as well
    private static double lengthOption2 = 2.953f;
    private static double lengthOption3 = 6.437f;

    private int N = 0; //used as top-level dimension of square matrix that the determinant is being calculated for

    private struct ke_kpe_Te_Vals{
        public double[,] kelem;
        public double[,] k_Local_elem;
        public double[,] Telem;
    }

    private double[,] U;
    private double[,] Pc;
    private int UDOF;
 
    #endregion

    #region Public Variables
    //These all have to be public so they can be accessed by reference, either by the UI buttons or by other classes (like ARManager). 

    /*
        Points on structure:
        B------------C
        |            |
        |            |
        A            D
    */
    public struct calHatPc_Vals{
        public double[,] finalValues;
        public int size;
    }

    public string configCase;

    public bool APinned;
    public bool BPinned;
    public bool CPinned;
    public bool DPinned;

    public double ABLength = lengthOption1;
    public double BCLength = lengthOption1;
    public double CDLength = lengthOption1;

    public int NDOF; //consider restructuring as private with getters and setters, since these will not be changed after configure is called
    public int Ncon;
    public int[,] le; //Note: Uses [row,col] structure.

    public GameObject Canvas1;
    public GameObject Canvas2;
    public GameObject NextButton;
    public GameObject ErrorText;

    public TMP_Dropdown ABDrop;
    public TMP_Dropdown BCDrop;
    public TMP_Dropdown CDDrop;

    #endregion

    // Start is called before the first frame update
    void Start()
    {
        DontDestroyOnLoad(this);

        //Using case 1 as default
        NDOF = 14;
        Ncon = 5;
        le = new int[3,6] {{10, 11, 9, 1, 2, 8}, {12, 13, 14, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
    }
    
    public void ToggleA(){
        APinned = !APinned;
    }

    public void ToggleB(){
        BPinned = !BPinned;
    }

    public void ToggleC(){
        CPinned = !CPinned;
    }

    public void ToggleD(){
        DPinned = !DPinned;
    }

    public void UpdateABLength(){
        int option = ABDrop.value;
        if(option == 0){
            ABLength = lengthOption1;

            CDLength = lengthOption1; //AB and CD need to be the same length

            CDDrop.value = ABDrop.value;

        } else if(option == 1){
            ABLength = lengthOption2;

            CDLength = lengthOption2; //AB and CD need to be the same length

            CDDrop.value = ABDrop.value;

        }else if(option == 2){
            ABLength = lengthOption3;

            CDLength = lengthOption3; //AB and CD need to be the same length

            CDDrop.value = ABDrop.value;
        }
    }

    public void UpdateBCLength(){
        int option = BCDrop.value;
        if(option == 0){
            BCLength = lengthOption1;

        } else if(option == 1){
            BCLength = lengthOption2;

        }else if(option == 2){
            BCLength = lengthOption3;
        }
    }

    public void UpdateCDLength(){
        int option = CDDrop.value;
        if(option == 0){
            ABLength = lengthOption1;

            CDLength = lengthOption1; //AB and CD need to be the same length

            ABDrop.value = CDDrop.value;


        } else if(option == 1){
            ABLength = lengthOption2;

            CDLength = lengthOption2; //AB and CD need to be the same length

            ABDrop.value = CDDrop.value;

        }else if(option == 2){
            ABLength = lengthOption3;
            
            CDLength = lengthOption3; //AB and CD need to be the same length

            ABDrop.value = CDDrop.value;
        }
    }

    public void CheckHowManyPinned(){
        if(APinned && BPinned && CPinned && DPinned){
            ErrorText.SetActive(true);
            NextButton.SetActive(false);
        } else {
            ErrorText.SetActive(false);
            NextButton.SetActive(true);
        }
    }

    public void UINext(){
        Canvas1.SetActive(false);
        Canvas2.SetActive(true);
    }

    public void UIBack(){
        Canvas1.SetActive(true);
        Canvas2.SetActive(false);
    }

    public void LoadMainScene(){
        Configure();
        U = CalReactions();
        calHatPc_Vals answer = calHatPc(1f/0.0254f);
        string answerStr = "";
        for(int i = 0; i < NDOF-UDOF; i++){
            answerStr = answerStr + "\n" + answer.finalValues[i,0];
        }        
        Debug.Log("SOLUTION: \n" + answerStr);
        SceneManager.LoadScene("MainARScene");
    }

    public void LoadGroundPlaneScene(){
        Configure();
        U = CalReactions();
        calHatPc_Vals answer = calHatPc(1f/0.0254f);
        string answerStr = "";
        for(int i = 0; i < NDOF-UDOF; i++){
            answerStr = answerStr + "\n" + answer.finalValues[i,0];
        }        
        Debug.Log("SOLUTION: \n" + answerStr);

        SceneManager.LoadScene("ARScenePlane");
    }

    public string calculateSolution(double displacement){
        calHatPc_Vals answer = calHatPc(displacement);
        string answerStr = "";
        for(int i = 0; i < NDOF-UDOF; i++){
            answerStr = answerStr + "\n" + answer.finalValues[i,0];
        }        
        Debug.Log("SOLUTION: \n" + answerStr);
        return answerStr;
    }

    public void LoadMenuScene(){
        SceneManager.LoadScene("MainMenu");
    }

    /*
        Points on structure:
        B------------C
        |            |
        |            |
        A            D
        p = pinned
        f = fixed
        Order of cases: ABCD
        EX: "pfpf" means A is pinned, B is fixed, C is pinned, D is fixed
    */
    void Configure(){
        string a = "f";
        string b = "f";
        string c = "f";
        string d = "f";
        if(APinned){ a = "p";}
        if(BPinned){ b = "p";}
        if(CPinned){ c = "p";}
        if(DPinned){ d = "p";}
        configCase = a + b + c + d;
        
        switch (configCase){
            case "pppf": // 1
                NDOF = 14;
                Ncon = 5;
                //le = new int[3,6] {{10, 11, 9, 1, 2, 8}, {12, 13, 14, 4, 5, 7}, {1, 2, 3, 4, 5, 6}}; //these should be transposed
                le = new int[6,3] { {9,11,0}, 
                                    {10,12,1}, 
                                    {8,13,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {7,6,5} };
                Debug.Log("CASE: 1"); 
                break;
            case "pfpp": // 2
                NDOF = 13;
                Ncon = 4;
                //le = new int[3,6] {{10, 11, 9, 1, 2, 3}, {12, 13, 8, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {9,11,0}, 
                                    {10,12,1},
                                    {8,7,2},
                                    {0,3,3}, 
                                    {1,4,5}, 
                                    {2,6,5} };
                Debug.Log("CASE: 2"); 
                break;
            case "ppfp": // 3
                NDOF = 13;
                Ncon = 4;
                // le = new int[3,6] {{10, 11, 9, 1, 2, 7}, {12, 13, 8, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {9,11,0}, 
                                    {10,12,1}, 
                                    {8,7,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {6,5,5} };
                Debug.Log("CASE: 3"); 
                break;            
            case "fppp": // 4
                NDOF = 14;
                Ncon = 5;
                // le = new int[3,6] {{10, 11, 9, 1, 2, 8}, {12, 13, 14, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {9,11,0}, 
                                    {10,12,1}, 
                                    {8,13,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {7,6,5} };
                Debug.Log("CASE: 4"); 
                break;            
            case "ffpp": // 5
                NDOF = 13;
                Ncon = 5;
                // le = new int[3,6] {{9, 10, 11, 1, 2, 3}, {12, 13, 8, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,11,0}, 
                                    {9,12,1}, 
                                    {10,7,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,6,5} };
                Debug.Log("CASE: 5"); 
                break;            
            case "fpfp": // 6
                NDOF = 13;
                Ncon = 5;
                // le = new int[3,6] {{9, 10, 11, 1, 2, 7}, {12, 13, 8, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,11,0}, 
                                    {9,12,1}, 
                                    {10,7,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {6,5,5} };
                Debug.Log("CASE: 6"); 
                break;            
            case "fppf": // 7
                NDOF = 14;
                Ncon = 6;
                // le = new int[3,6] {{9, 10, 11, 1, 2, 8}, {12, 13, 14, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,11,0}, 
                                    {9,12,1}, 
                                    {10,13,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {7,6,5} };
                Debug.Log("CASE: 7"); 
                break;            
            case "ppff": // 8
                NDOF = 13;
                Ncon = 5;
                // le = new int[3,6] {{9, 10, 8, 1, 2, 7}, {11, 12, 13, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,10,0}, 
                                    {9,11,1}, 
                                    {7,12,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {6,5,5} };
                Debug.Log("CASE: 8"); 
                break;
            case "fffp": // 9
                NDOF = 12;
                Ncon = 5;
                // le = new int[3,6] {{8, 9, 10, 1, 2, 3}, {11, 12, 7, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {7,10,0}, 
                                    {8,11,1}, 
                                    {9,6,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,5,5} };
                Debug.Log("CASE: 9"); 
                break;
            case "pfff": // 10
                NDOF = 12;
                Ncon = 5;
                // le = new int[3,6] {{8, 9, 7, 1, 2, 3}, {10, 11, 12, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {7,9,0}, 
                                    {8,10,1}, 
                                    {6,11,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,5,5} };
                Debug.Log("CASE: 10"); 
                break;
            case "fpff": // 11
                NDOF = 13;
                Ncon = 6;
                // le = new int[3,6] {{8, 9, 10, 1, 2, 7}, {11, 12, 13, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {7,10,0}, 
                                    {8,11,1}, 
                                    {9,12,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {6,5,5} };
                Debug.Log("CASE: 11"); 
                break;
            case "ffpf": // 12
                NDOF = 13;
                Ncon = 6;
                // le = new int[3,6] {{8, 9, 10, 1, 2, 3}, {11, 12, 13, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {7,10,0}, 
                                    {8,11,1}, 
                                    {9,12,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,6,5} };
                Debug.Log("CASE: 12"); 
                break;
            case "ffff": // 13
                NDOF = 12;
                Ncon = 6;
                // le = new int[3,6] {{7, 8, 9, 1, 2, 3}, {10, 11, 12, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {6,9,0}, 
                                    {7,10,1}, 
                                    {8,11,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,5,5} };
                Debug.Log("CASE: 13"); 
                break;
            case "pfpf": // 14
                NDOF = 13;
                Ncon = 5;
                // le = new int[3,6] {{9, 10, 8, 1, 2, 3}, {11, 12, 13, 4, 5, 7}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,10,0}, 
                                    {9,11,1}, 
                                    {7,12,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,6,5} };
                Debug.Log("CASE: 14"); 
                break;
            case "pffp": // 15
                NDOF = 12;
                Ncon = 4;
                // le = new int[3,6] {{9, 10, 8, 1, 2, 3}, {11, 12, 7, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
                le = new int[6,3] { {8,10,0}, 
                                    {9,11,1}, 
                                    {7,6,2}, 
                                    {0,3,3}, 
                                    {1,4,4}, 
                                    {2,5,5} };
                Debug.Log("CASE: 15"); 
                break;
            case "pppp": // 16
                Debug.Log("ERROR: Case 16 is not allowed");
                break;
            default:
                Debug.Log("ERROR: Invalid case");
                break;
        }
    }

    double[,] CalReactions()
    {
        double[] lengths = { ABLength, BCLength, CDLength };
        int Nelem = 3;
        int EDOF = 6;
        double[] Angle = { 90 * Mathf.PI / 180, 90 * Mathf.PI / 180, 0 };
        double[] AE = { 10.1, 10.1, 10.1 };
        double[] EI = { 0.711442913, 0.711442913, 0.711442913 };

        //P = [1 zeros(1,NDOF-1)]';     % one unit load
        //float[] P = {1};
        List<double> P = new List<double>();
        P.Add(1);
        for (int i = 1; i < NDOF-1; i++)
        {
            P.Add(0);
        }

        //Uc = zeros(1,Ncon)';          % in.
        //Uc would equal [0,0,0,0,0...], but this is asking for the transpose of that
        // List<double> Uc = new List<double>();
        // for(int i = 0; i < Ncon; i++)
        // {
        //     Uc.Add(0);
        // }

        double[,] Uc = new double[Ncon,1];
        for(int i = 0; i < Ncon; i++)
        {
            Uc[i,0] = 0;
        }

        //K = zeros(NDOF, NDOF);
        //double[,] K = {};
        double[,] K = new double[NDOF,NDOF];
        for(int i = 0; i < NDOF; i++){
            for(int j = 0; j < NDOF; j++){
                K[i,j] = 0;
            }
        }

        /*
         for e = 1:Nelem
            [ke1,kpe1, Te1] = k_frame_elem2(EI(e),Length(e),AE(e),Angle(e));
            K11=K_SMADM(ke1,le(:,e),NDOF); %Note, le(:,e) means the eth column of le
            K = K + K11;
         end
         */
        for (int e = 0; e < Nelem; e++)
        {
            ke_kpe_Te_Vals ke_kpe_Te = k_frame_elem2(EI[e], lengths[e], AE[e], Angle[e]);

            //Get the eth column of le
            int[] leCol = new int[6];
            for(int i = 0; i < EDOF; i++){
               leCol[i] = le[i,e];
            }

            double[,] K11 = K_SMADM(ke_kpe_Te.kelem, leCol);

            Debug.Log("PRINTING K11:");
            printSquareMatrix(K11,NDOF);

            //K = K + K11;
            for(int i = 0; i < NDOF; i++){
                for(int j = 0; j < NDOF; j++){
                    K[i,j] += K11[i,j];
                    // Debug.Log(K[i][j] + " ADDED TO K["+i+"]["+j+"]");
                }
            }
        }

        //Every item here is a matrix that is partitioned from K:
        
        //UDOF = NDOF-Ncon;
        UDOF = NDOF - Ncon;

       //Kuu  = K(1:UDOF,1:UDOF);
        double[,] Kuu = new double[UDOF,UDOF];
        for(int i = 0; i < UDOF; i++){
            for(int j = 0; j < UDOF; j++){
                Kuu[i,j] = K[i,j];
            }
        }
        Debug.Log("KUU MATRIX:\n");
        printSquareMatrix(Kuu, UDOF);
        
        //Kcu  = K(UDOF+1:NDOF,1:UDOF);
        double[,] Kcu = new double[NDOF-UDOF, UDOF];
        for(int i = 0; i < NDOF-UDOF; i++){
            for(int j = 0; j < UDOF; j++){
                Kcu[i,j] = K[UDOF+i,j];
            }
        }
        
        //Kuc  = K(1:UDOF,UDOF+1:NDOF);
        double[,] Kuc = new double[UDOF, NDOF-UDOF];
        for(int i = 0; i < UDOF; i++){
            for(int j = 0; j < NDOF-UDOF; j++){
                Kuc[i,j] = K[i,UDOF+j];
            }
        }
        
        //Kcc  = K(UDOF+1:NDOF,UDOF+1:NDOF);
        double[,] Kcc = new double[NDOF-UDOF, NDOF-UDOF];
        for(int i = 0; i < NDOF-UDOF; i++){
            for(int j = 0; j < NDOF-UDOF; j++){
                Kcu[i,j] = K[UDOF+i,UDOF+j];
            }
        }

        //P = [1 zeros(1,NDOF-1)]';     % one unit load, note that transpose is taken to make it one column instead of one row
        //Pu   = P(1:UDOF,1);
       double[,] Pu = new double[UDOF,1];
        for(int i = 0; i < UDOF; i++){
            Pu[i,0] = P[i];
        }

        //Uu = Kuu\(Pu - Kuc * Uc);
        // //In order to get the inverse of Kuu, we need its determinant
        // //SET N here as size of square matrix kuu;
        N = UDOF;
        double det = determinantOfMatrix(Kuu, N);
        Debug.Log("DETERMINANT: " + det);

        //Now find the inverse of Kuu
        double[,] invKuu = new double[UDOF, UDOF];
        if(!inverse(Kuu, invKuu,N)){
            Debug.Log("FAIL: Cannot find inverse of matrix Kuu!");
            double[,] falseReturn = {{-9999,-9999,-9999}};
            return falseReturn;
        } else{
            Debug.Log("SUCCESS: Inverse of matrix Kuu found!");
        }

        //Kuc is UDOF x NDOF-UDOF
        //Uc is a list of size Ncon, so 1xNcon, transposed, so Ncon x 1
        //Pu is a list of size UDOF, matrix is UDOF x 1
        //Example - Case 1:
        //    NDOF = 14;
        //    Ncon = 5;
        //    UDOF = NDOF - NCON = 9; 
        //    double[,] Kuc = new double[UDOF, NDOF-UDOF];
        //    Kuc = 9 x 5
        //    Uc = 5 x 1
        //    Result of multiplication: 9 x 1
        //    Pu - (Kuc*Uc) = (9x1)-(9x1)
        //    Kuu\ or "invKuu" is UDOF x UDOF so 9 x 9
        //    (9 x 9) * (9 x 1) = 9 x 1

        double[,] KucUcProduct = multiplyMatrices(Kuc, Uc, UDOF, NDOF-UDOF, Ncon, 1);
        double[,] PuMinusProduct = new double[UDOF, 1];
        for(int i = 0; i < UDOF; i++){
            PuMinusProduct[i,0] = Pu[i,0] - KucUcProduct[i,0]; 
        }
        
        double[,] Uu = multiplyMatrices(invKuu, PuMinusProduct, UDOF, UDOF, UDOF, 1);
        
        // Pc = Kcu*Uu + Kcc*Uc;               % Uu and Uc are vectors, Kcu and Kcc are matrices
        //Uu is 9 x 1
        //Kcu is NDOF-UDOF x UDOF = 5 x 14
        double[,] KcuUUProduct = multiplyMatrices(Kcu, Uu, NDOF-UDOF, UDOF, UDOF, 1); //result = NDOF-UDOF x 1
        double[,] KccUcProduct = multiplyMatrices(Kcc, Uc, NDOF-UDOF, NDOF-UDOF, Ncon,1); //result = NDOF-UDOF x 1
        Pc = new double[NDOF-UDOF, 1];
        for(int i = 0; i < NDOF-UDOF; i++){
            Pc[i,0] = KcuUUProduct[i,0] + KccUcProduct[i,0];
        }

        // U(1:UDOF,1)=Uu;                     % U doesn't exist until here
        // U(UDOF+1:NDOF,1)=Uc;                % this part of the matrix isn't used
        double[,] U = new double[NDOF,1];
        for(int i = 0; i < UDOF; i++){
            U[i,0] = Uu[i,0];
        }
        for(int i = 0; i < NDOF-UDOF; i++){
            U[i+UDOF,0] = Uc[i,0];
        }
        return U;
    }

    public calHatPc_Vals calHatPc(double displacement){

        //NOTE: Displacement is received in m/cm, must convert to inches
        displacement *= 0.0254f; //conversion from meters to inches

        /*
        haU1 = U(1, 1);                        % Displacement measurement at the specific location(this function need be deleted)
        hatPc = (haU1 / U(1, 1)) * Pc          % hatPc is a vector, one single column of the matrix but the number of items may change, maybe store as matrix
        */

        //Pc = new double[NDOF-UDOF, 1];

        double[,] hatPc = new double[NDOF-UDOF,1];
        for(int i = 0; i < NDOF-UDOF; i++){
            //hatPc[i,0] = (displacement/U[0,0]) * Pc[i,0]; //use this when displacement is achieved
            hatPc[i,0] = (displacement/U[0,0]) * Pc[i,0];
        }
        calHatPc_Vals finalReturn;
        finalReturn.finalValues = hatPc;
        finalReturn.size = NDOF-UDOF;
        return finalReturn;
    }

    /*
    function [kelem, k_Local_elem, Telem]=k_frame_elem2(EI,L,AE,Theta)
    
    AE  =  AE/L;                   % these are numbers, store as double/float?
    EI12 = 12*EI/L^3;
    EI6 = 6*EI/L^2;
    EI4 = 4*EI/L;
    EI2 = 2*EI/L;
    k_Local_elem  = [   AE    0    0   -AE    0    0     ;
                        0  EI12   EI6    0 -EI12   EI6   ;
                        0   EI6   EI4    0  -EI6   EI2   ;
                        -AE    0    0    AE    0    0     ;
                        0 -EI12  -EI6    0  EI12  -EI6   ;
                        0   EI6   EI2    0  -EI6   EI4  ];
    sine = sin(Theta);  cosine = cos(Theta);
    Telem = [  cosine    sine    0     0       0    0   ;
                -sine    cosine   0     0       0    0   ;
                    0        0      1     0       0    0   ;
                    0        0      0   cosine  sine   0   ;
                    0        0      0  -sine  cosine   0   ;
                    0        0      0    0        0    1  ];
    kelem = Telem' * k_Local_elem * Telem;
    */
    ke_kpe_Te_Vals k_frame_elem2(double EI, double Length , double AE, double Angle )
    {
        double[,] kelem = new double[6,6];

        AE = AE/Length;
        double EI12 = (12 * EI)/(Length * Length * Length);
        double EI6 = (6*EI)/(Length * Length);
        double EI4 = (4 * EI)/Length;
        double EI2 = (2 * EI)/Length;
        //Debug.Log("EI12: " + EI12 + "\nEI6" + EI6 + "\nEI4" + EI4 + "\nEI2" + EI2);       //Confirmed these have nonzero values 

        double[,] k_Local_elem = 
                            {{AE, 0, 0, AE * -1, 0, 0}, 
                            {0, EI12, EI6, 0, EI12 * -1, EI6},
                            {0, EI6, EI4, 0, EI6 * -1, EI2}, 
                            {AE * -1, 0, 0, AE, 0, 0}, 
                            {0, EI12 * -1, EI6 * -1, 0, EI12, EI6 * -1}, 
                            {0, EI6, EI2, 0, EI6 * -1, EI4} };

        double sine = Mathf.Sin((float)Angle);
        double cosine = Mathf.Cos((float)Angle);
        
        double[,] Telem = 
                            {{cosine, sine, 0, 0, 0, 0}, 
                            {sine * -1, cosine, 0, 0, 0, 0},
                            {0, 0, 1, 0, 0, 0}, 
                            {0, 0, 0, cosine, sine, 0}, 
                            {0, 0, 0, sine * -1, cosine, 0}, 
                            {0, 0, 0, 0, 0, 1}};

        double[,] TelemTranspose = 
                            {{cosine, sine * -1, 0, 0, 0, 0}, 
                            {sine, cosine, 0, 0, 0, 0},
                            {0, 0, 1, 0, 0, 0}, 
                            {0, 0, 0, cosine, sine * -1, 0}, 
                            {0, 0, 0, sine, cosine, 0}, 
                            {0, 0, 0, 0, 0, 1}};

        kelem = multiplyMatrices(TelemTranspose, k_Local_elem, 6, 6, 6, 6);
        kelem = multiplyMatrices(kelem, Telem, 6,6,6,6);

        ke_kpe_Te_Vals toReturn;
        toReturn.kelem = kelem;
        toReturn.k_Local_elem = k_Local_elem;
        toReturn.Telem = Telem;
        Debug.Log("k_frame_elem2 RETURN VALUES:");
        Debug.Log("kelem:");
        printSquareMatrix(kelem, 6);
        Debug.Log("k_Local_elem");
        printSquareMatrix(k_Local_elem,6);
        Debug.Log("Telem:");
        printSquareMatrix(Telem,6);
        return toReturn;
    }

    /*
    K=K_SMADM(ke,le,DOF)
    % K Stiffness Matix Assembly using Direct Method  SMADM
    % Sdiq Anwar Taher and Dr. Li December 27, 2019
    % The University of Kansas
    K = zeros(DOF,DOF);
    edof=size(le,1);
    %
    for ii = 1:edof
        i = le(ii);
        if i > 0
            for jj = 1:edof
                j = le(jj);
                if j > 0
                    K(i,j) = ke(ii,jj);
                end
            end
        end
    end
    */

    double[,] K_SMADM(double[,] ke1, int[] le)
    {
        double[,] K = new double[NDOF,NDOF];
        for(int i = 0; i < NDOF; i++){
            for(int j = 0; j < NDOF; j++){
                K[i,j] = 0;
            }
        }

        int edof = le.Length;

        for(int i = 0; i < edof; i++){
            int leI = le[i];
            if(leI >= 0){
                for(int j = 0; j < edof; j++){
                    int leJ = le[j];
                    if(leJ >= 0){
                        K[leI,leJ] = ke1[i,j];
                    }
                }
            }
        }

        return K;
    }

    //***Below functions are taken from GeeksForGeeks: https://www.geeksforgeeks.org/determinant-of-a-matrix/***
      
    // Function to get cofactor of  
    // mat[p][q] in temp[][]. n is  
    // current dimension of mat[][] 
    void getCofactor(double [,]mat, double [,]temp, int p, int q, int n) 
    { 
        int i = 0, j = 0; 
      
        // Looping for each element of  
        // the matrix 
        for (int row = 0; row < n; row++) 
        { 
            for (int col = 0; col < n; col++) 
            { 
                  
                // Copying into temporary matrix  
                // only those element which are  
                // not in given row and column 
                if (row != p && col != q) 
                { 
                    temp[i, j++] = mat[row, col]; 
      
                    // Row is filled, so increase  
                    // row index and reset col  
                    //index 
                    if (j == n - 1) 
                    { 
                        j = 0; 
                        i++; 
                    } 
                } 
            } 
        } 
    } 
      
    /* Recursive function for 
       finding determinant 
       of matrix. n is current  
       dimension of mat[][]. */
    double determinantOfMatrix(double [,]mat, int n) 
    { 
        double D = 0; // Initialize result 
      
        // Base case : if matrix  
        // contains single 
        // element 
        if (n == 1) 
            return mat[0, 0]; 
          
        // To store cofactors 
        double [,]temp = new double[N, N];  
          
        // To store sign multiplier 
        int sign = 1;  
      
        // Iterate for each element 
        // of first row 
        for (int f = 0; f < n; f++) 
        { 
              
            // Getting Cofactor of mat[0][f] 
            getCofactor(mat, temp, 0, f, n); 
            D += sign * mat[0, f]  
            * determinantOfMatrix(temp, n - 1); 
      
            // terms are to be added with  
            // alternate sign 
            sign = -sign; 
        } 
      
        return D; 
    } 

    // Function to get adjoint of A[N,N] in adj[N,N]. 
    void adjoint(double [,]A, double [,]adj, int N) 
    { 
        if (N == 1) 
        { 
            adj[0, 0] = 1; 
            return; 
        } 
    
        // temp is used to store cofactors of [,]A 
        int sign = 1; 
        double [,]temp = new double[N, N]; 
    
        for (int i = 0; i < N; i++) 
        { 
            for (int j = 0; j < N; j++) 
            { 
                // Get cofactor of A[i,j] 
                getCofactor(A, temp, i, j, N); 
    
                // sign of adj[j,i] positive if sum of row 
                // and column indexes is even. 
                sign = ((i + j) % 2 == 0)? 1: -1; 
    
                // Interchanging rows and columns to get the 
                // transpose of the cofactor matrix 
                adj[j, i] = (sign) * (determinantOfMatrix(temp, N - 1)); 
            } 
        } 
    } 
  
    // Function to calculate and store inverse, returns false if 
    // matrix is singular 
    bool inverse(double [,]A, double [,]inverse, int N) 
    { 
        // Find determinant of [,]A 
        double det = determinantOfMatrix(A, N); 
        if (det == 0) 
        { 
            Debug.Log("Singular matrix, can't find its inverse"); 
            return false; 
        } 
    
        // Find adjoint 
        double [,]adj = new double[N, N]; 
        adjoint(A, adj, N); 
    
        // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
        for (int i = 0; i < N; i++) 
            for (int j = 0; j < N; j++) 
                inverse[i, j] = adj[i, j]/det; 
    
        return true; 
    }

    void printSquareMatrix(double[,] M, int N){
        string printString = "";
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                printString = printString + M[i,j] + ", ";
            }
            printString += "\n";
        }
        Debug.Log(printString);
    } 

    //Below function adapted from tutorialspoint: https://www.tutorialspoint.com/chash-program-to-multiply-two-matrices

    double[,] multiplyMatrices(double[,] M1, double[,] M2, int m, int n, int p, int q){

        double[,] mult = new double[m, q];

         if(n != p) {
            Debug.Log("ERROR: Matrix multiplication not possible");
         } else {
            for (int i = 0; i < m; i++) {
               for (int j = 0; j < q; j++) {
                  mult[i, j] = 0;
                  for (int k = 0; k < n; k++) {
                     mult[i, j] += M1[i, k] * M2[k, j];
                  }
               }
            }
        }
        return mult;
    }
}
