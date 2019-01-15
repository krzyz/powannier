#ifndef POWANNIER_WANNIERPOSITIONS_H
#define POWANNIER_WANNIERPOSITIONS_H

#include <vector>

namespace POWannier {
  /**
   * @brief Specification of the degeneracy of the eigenvalues of
   * \f$ \hat{\lambda_i} \f$ operator.
   * 
   * While diagonalizing the operators \f$ \hat{\lambda_i} \f$ (for their
   * definitions see RSystem class) some of the eigenvalues may be degenerate
   * in the single elementary cell. Because of that, when choosing the subspace
   * of the degenerated eigenvalues (second operation in each diagonalization
   * loop step), this degeneracy has to be included. The structure of 
   * WannierPositions is <tt>k</tt>-ary tree (with @p k equal to the number
   * of bands considered, which is also equal to the number of %Wannier
   * functions constructed from these bands in each elementary cell) of 
   * the height @p dim (dimension of the problem). The number of <tt>i</tt>-th
   * level children of some <tt>i-1</tt>-th level node corresponds to the
   * number of distinct eigenvalues (modulo 1, as this is the width of the\
   * elementary cell in units of \f$ \hat{\lambda_i} \f$) when diagonalizing 
   * the \f$ \hat{\lambda}_i \f$ operator in the subspace chosen by 
   * @p wanLocation corresponding to that parent mode.
   * 
   * In order to better understand how to structure the tree, using one of the
   * constructors WannierPositions(), a few examples follow. Each of the cases
   * assumes some spatial composition of %Wannier functions which is represented
   * in the middle column. In the left column few possible and equivalent 
   * choices of arguments of WannierPositions() are shown, each recreating 
   * the tree structure shown in the right column (the leaves of which 
   * correspond to the positions of %Wannier functions, the set of indices 
   * which which are listed there are the appropriate values of @p wanLocation
   * (see RSystem::getWannier()).
   * The following simplifications can be made to the arguments provided
   * to WannierPositions():
   * - \f$ \underbrace{\{1,1,\dots, 1\}}_{n} = \{ n \} \f$
   * - \f$ \{\{\dots \{1\} \dots \} \} = 1 \f$.
   * 
   * @par 2D example
   * <table style="{border:none;}"><tr><td>
   * <tt>{{1},{1,1},{1}}</tt> \n
   * <tt>{1,2,1} </tt>
   * </td><td>
   * @image html 2d.svg
   * </td><td>
   * @image html tree1.svg
   * </td></tr></table>
   * @par 3D example
   * <table style="{border:none;}"><tr><td>
   * <tt>{{{1,1},{1,1}},{{1}},{{1,1},{1,1}}</tt> \n
   * <tt>{{2,2},{1},{2,2}}</tt> \n
   * <tt>{{2,2},1,{2,2}}</tt>
   * </td><td>
   * @image html 3d.svg
   * </td><td>
   * @image html tree2.svg
   * </td></tr></table>

   * @see RSystem
   */
  class WannierPositions {
    public:
      /// Default constructor;
      WannierPositions() {};

      /**
       * @brief Constructor using initializer list of ints.
       * @param leafNs
       *   List of nodes, where each node contains only leaves and 
       *   is represented by an integer number.
       */
      WannierPositions(
          std::initializer_list<int> leafNs);

      /**
       * @brief Constructor using @p std::vector of ints.
       * @param leafNs
       *   List of nodes, where each node contains only leaves and 
       *   is represented by an integer number.
       */
      WannierPositions(
          std::vector<int> leafNs);

      /**
       * @brief Constractor using a list of children nodes.
       * @param children
       *   List of WannierPosition nodes. Because of implicit conversion
       *   of initializer lists (WannierPositions(std::initializer_list<int>))
       *   this constructor can recursively add nodes specified by \{...\} 
       *   lists, such as {{2,2},1,{2,2}}.
       */
      WannierPositions(
          std::initializer_list<WannierPositions> children);

      /**
       * @brief Add another node as a child.
       * @param child
       *   A child node to be added.
       */
      void addChild(WannierPositions child);

      /**
       * @brief Get the number of children.
       * @returns
       *   The number of children.
       */
      int childrenNumber() const;

      /**
       * @brief Get the number of descendants that are leaves.
       * @returns
       *   The number of descending leaves.
       */
      int descendantLeavesNumber() const;

      /**
       * @brief Get sum of leaves descending from all childrens that are
       *   on the left of a given child.
       * @param i
       *   Index of the child.
       * @returns
       *   The number of leaves descending from the children that are
       *   on the left of the child indexed by <tt>i</tt> (starts with 0).
       */
      int descendantLeavesNumberLeftTo(int i) const;

      /**
       * @brief Get the <tt>i</tt>-th child node.
       * @param i The index of the children node to be return.
       * @returns The <tt>i</tt>-th child node.
       */
      WannierPositions getChild(int i);

    private:
      template <class It>
      WannierPositions(It first, It last);

      std::vector<WannierPositions> _children;
  };
}

#endif