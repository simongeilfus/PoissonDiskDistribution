/*
 Copyright (c) 2015 Simon Geilfus
 
 Algorithm from Fast Poisson Disk Sampling in Arbitrary Dimensions by Robert Bridson
 http://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
 as explained in this article: http://devmag.org.za/2009/05/03/poisson-disk-sampling/
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#include "PoissonDiskDistribution.h"

#include "cinder/Rand.h"
#include "cinder/Log.h"

using namespace std;
using namespace ci;

namespace {	
	class Grid {
	public:
		using Cell = std::vector<vec2>;
		
		Grid( const ci::Area &bounds, uint32_t k );
		
		void add( const ci::vec2 &position );
		bool hasNeighbors( const ci::vec2 &p, float radius );
		
		void resize( const ci::Area &bounds, uint32_t k );
		void resize( uint32_t k );
		
	protected:
		std::vector<Cell>	mGrid;
		ci::ivec2			mNumCells, mOffset;
		ci::Area			mBounds;
		uint32_t			mK, mCellSize;
	};
	
	
	Grid::Grid( const Area &bounds, uint32_t k )
	{
		resize( bounds, k );
	}
	
	void Grid::add( const ci::vec2 &position )
	{
		int x = ( ( uint32_t ) position.x + mOffset.x ) >> mK;
		int y = ( ( uint32_t ) position.y + mOffset.y ) >> mK;
		int j = x + mNumCells.x*y;
		
		if( j < mGrid.size() ){
			mGrid[j].push_back( position );
		}
		else {
			CI_LOG_E( "Out of bounds" );
		}
	}
	
	bool Grid::hasNeighbors( const ci::vec2 &p, float radius )
	{
		float sqRadius	= radius * radius;
		ivec2 radiusVec = ivec2( radius );
		ivec2 min       = glm::max( glm::min( ivec2( p ) - radiusVec, ivec2( mBounds.getLR() ) - ivec2( 1 ) ), ivec2( mBounds.getUL() ) );
		ivec2 max       = glm::max( glm::min( ivec2( p ) + radiusVec, ivec2( mBounds.getLR() ) - ivec2( 1 ) ), ivec2( mBounds.getUL() ) );
		
		ivec2 minCell( ( min.x + mOffset.x ) >> mK, ( min.y + mOffset.y ) >> mK );
		ivec2 maxCell = glm::min( 1 + ivec2( ( max.x + mOffset.x ) >> mK, ( max.y + mOffset.y ) >> mK ), mNumCells );
		for( unsigned y = minCell.y; y < maxCell.y; y++ ) {
			for( unsigned x = minCell.x; x < maxCell.x; x++ ) {
				for( auto cell : mGrid[ x + mNumCells.x*y ] ) {
					if( glm::length2( p - cell ) < sqRadius ){
						return true;
					}
				}
			}
		}
		return false;
	}
	
	void Grid::resize( const Area &bounds, uint32_t k )
	{
		mBounds = bounds;
		resize( k );
	}
	void Grid::resize( uint32_t k )
	{
		mK          = k;
		mCellSize   = 1 << k;
		mOffset     = glm::abs( mBounds.getUL() );
		mNumCells   = glm::ceil( vec2( mBounds.getSize() ) / (float) mCellSize );
		mGrid.clear();
		mGrid.resize( mNumCells.x * mNumCells.y );
	}
} // anonymous namespace

std::vector<glm::vec2> poissonDiskDistribution( float separation, const ci::Rectf &bounds, const std::vector<glm::vec2> &initialSet, int k )
{
	// prepare working structures
	vector<vec2> processingList;
	vector<vec2> outputList;
	
	// create grid
	Grid grid( Area( bounds ), 3 );
	
	// add the initial points
	for( auto p : initialSet ){
		processingList.push_back( p );
		outputList.push_back( p );
		grid.add( p );
	}
	
	// if there's no initial points add the center point
	if( !processingList.size() ){
		processingList.push_back( bounds.getCenter() );
		outputList.push_back( bounds.getCenter() );
		grid.add( bounds.getCenter() );
	}
	
	// while there's points in the processing list
	while( processingList.size() ){
		
		// pick a random point in the processing list
		int randPoint = randInt( 0, processingList.size() - 1 );
		vec2 center = processingList[randPoint];
		
		// remove it
		processingList.erase( processingList.begin() + randPoint );
		
		// spawn k points in an anulus around that point
		// the higher k is, the higher the packing will be and slower the algorithm
		for( int i = 0; i < k; i++ ){
			float randRadius	= separation * ( 1.0f + randFloat() );
			float randAngle		= randFloat() * M_PI * 2.0f;
			vec2 newPoint		= center + vec2( cos( randAngle ), sin( randAngle ) ) * randRadius;
			
			// check if the new random point is in the window bounds
			// and if it has no neighbors that are too close to them
			if( bounds.contains( newPoint )
			   && !grid.hasNeighbors( newPoint, separation ) ){
				
				// if the point has no close neighbors add it to the processing list, output list and grid
				processingList.push_back( newPoint );
				outputList.push_back( newPoint );
				grid.add( newPoint );
			}
		}
	}
	
	return outputList;
}

std::vector<glm::vec2> poissonDiskDistribution( const std::function<float(const glm::vec2&)> &distFunction, const ci::Rectf &bounds, const std::vector<glm::vec2> &initialSet, int k )
{
	// prepare working structures
	vector<vec2> processingList;
	vector<vec2> outputList;
	
	// create grid
	Grid grid( Area( bounds ), 3 );
	
	// add the initial points
	for( auto p : initialSet ){
		processingList.push_back( p );
		outputList.push_back( p );
		grid.add( p );
	}
	
	// if there's no initial points add the center point
	if( !processingList.size() ){
		processingList.push_back( bounds.getCenter() );
		outputList.push_back( bounds.getCenter() );
		grid.add( bounds.getCenter() );
	}
	
	// while there's points in the processing list
	while( processingList.size() ){
		
		// pick a random point in the processing list
		int randPoint = randInt( 0, processingList.size() - 1 );
		vec2 center = processingList[randPoint];
		
		// remove it
		processingList.erase( processingList.begin() + randPoint );
		
		// get the current min distance
		float dist = distFunction( center );
		
		// spawn k points in an anulus around that point
		// the higher k is, the higher the packing will be and slower the algorithm
		for( int i = 0; i < k; i++ ){
			float randRadius	= dist * ( 1.0f + randFloat() );
			float randAngle		= randFloat() * M_PI * 2.0f;
			vec2 newPoint		= center + vec2( cos( randAngle ), sin( randAngle ) ) * randRadius;
			
			// check if the new random point is in the window bounds
			// and if it has no neighbors that are too close to them
			if( bounds.contains( newPoint )
			   && !grid.hasNeighbors( newPoint, dist ) ){
				
				// if the point has no close neighbors add it to the processing list, output list and grid
				processingList.push_back( newPoint );
				outputList.push_back( newPoint );
				grid.add( newPoint );
			}
		}
	}
	
	return outputList;
}

std::vector<glm::vec2> poissonDiskDistribution( const std::function<float(const glm::vec2&)> &distFunction, const std::function<bool(const glm::vec2&)> &boundsFunction, const ci::Rectf &bounds, const std::vector<glm::vec2> &initialSet, int k )
{
	// prepare working structures
	vector<vec2> processingList;
	vector<vec2> outputList;
	
	// create grid
	Grid grid( Area( bounds ), 3 );
	
	// add the initial points
	for( auto p : initialSet ){
		processingList.push_back( p );
		outputList.push_back( p );
		grid.add( p );
	}
	
	// if there's no initial points add the center point
	if( !processingList.size() ){
		processingList.push_back( bounds.getCenter() );
		outputList.push_back( bounds.getCenter() );
		grid.add( bounds.getCenter() );
	}
	
	// while there's points in the processing list
	while( processingList.size() ){
		
		// pick a random point in the processing list
		int randPoint = randInt( 0, processingList.size() - 1 );
		vec2 center = processingList[randPoint];
		
		// remove it
		processingList.erase( processingList.begin() + randPoint );
		
		// get the current min distance
		float dist = distFunction( center );
		
		// spawn k points in an anulus around that point
		// the higher k is, the higher the packing will be and slower the algorithm
		for( int i = 0; i < k; i++ ){
			float randRadius	= dist * ( 1.0f + randFloat() );
			float randAngle		= randFloat() * M_PI * 2.0f;
			vec2 newPoint		= center + vec2( cos( randAngle ), sin( randAngle ) ) * randRadius;
			
			// check if the new random point is in the window bounds
			// and if it has no neighbors that are too close to them
			if( bounds.contains( newPoint )
			   && boundsFunction( newPoint )
			   && !grid.hasNeighbors( newPoint, dist ) ){
				
				// if the point has no close neighbors add it to the processing list, output list and grid
				processingList.push_back( newPoint );
				outputList.push_back( newPoint );
				grid.add( newPoint );
			}
		}
	}
	
	return outputList;
}