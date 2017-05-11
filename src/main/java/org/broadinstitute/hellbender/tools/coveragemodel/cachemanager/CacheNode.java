package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;

/**
 * The base class for all cache nodes in an {@link ImmutableComputableGraph}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
abstract class CacheNode {
    /**
     * A string identifier for the cache node
     */
    private final String key;

    /**
     * The collection of string identifiers of the immediate parents of this node (can be empty)
     */
    private final Collection<String> parents;

    /**
     * The collection of string identifiers of the tags associated to this node (can be empty)
     */
    private final Collection<String> tags;

    /**
     * Public constructor
     *
     * @param key string identifier of the cache node
     * @param tags the tags associated to this cache node
     * @param parents immediate parents of this cache node
     */
    CacheNode(@Nonnull final String key,
              @Nonnull final Collection<String> tags,
              @Nonnull final Collection<String> parents) {
        this.key = Utils.nonNull(key, "The key of a cache node can not be null");
        this.tags = Collections.unmodifiableCollection(Utils.nonNull(tags, "The tag collection of a cache node can not be null"));
        this.parents = Collections.unmodifiableCollection(Utils.nonNull(parents, "The immediate parents of a cache node can not be null"));
    }

    /**
     * Get the value stored in the node
     *
     * @param parents parent values (as a map from their string identifiers to their values)
     * @return a {@link Duplicable}; possibly by reference
     */
    abstract Duplicable get(@Nonnull final Map<String, Duplicable> parents);

    /**
     * Set the value of the node
     *
     * @param newValue new value; possibly stored by reference
     * @throws UnsupportedOperationException if the node is automatically computable
     */
    abstract void set(@Nullable final Duplicable newValue) throws UnsupportedOperationException;

    /**
     * Is the node primitive?
     */
    abstract boolean isPrimitive();

    /**
     * Is the node initialized yet?
     */
    abstract boolean hasValue();

    /**
     * Is the node externally computed?
     */
    abstract boolean isExternallyComputed();

    /**
     * Duplicate the node with updated value
     *
     * @param newValue new value; possibly stored by reference
     * @return a new {@link CacheNode} with the same key, parents, and tags but with a new value
     * @throws UnsupportedOperationException if the node is automatically computable
     */
    abstract CacheNode duplicateWithUpdatedValue(final Duplicable newValue) throws UnsupportedOperationException;

    /**
     * Make a deep copy of the node
     *
     * @return a deeply copied instance of {@link CacheNode}
     */
    abstract CacheNode duplicate();

    /**
     * Get the string identifier of the node
     * @return a non-null {@link String}
     */
    final String getKey() {
        return key;
    }

    /**
     * Get the collection of string identifier of the parents of this node (can be empty)
     */
    final Collection<String> getParents() {
        return Collections.unmodifiableCollection(parents);
    }

    /**
     * Get the collection of string identifier of the tags associated to this node (can be empty)
     */
    final Collection<String> getTags() {
        return Collections.unmodifiableCollection(tags);
    }

    @Override
    public final String toString() {
        return key;
    }

    /**
     * NOTE: equality comparison is done just based on the key
     * @param other another object
     */
    @Override
    public final boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;
        return (key.equals(((CacheNode) other).key));
    }

    /**
     * NOTE: hashcode is generated just based on the key
     */
    @Override
    public final int hashCode() {
        return key.hashCode();
    }
}
